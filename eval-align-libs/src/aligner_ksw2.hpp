/*
    SPEAR - Sorting Petabytes of Environmental and Ancient Reads
    Copyright (C) 2026 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#pragma once

#include "stats.hpp"

extern "C" {
#include "ksw2.h"
}

#include <cstdint>
#include <string>
#include <vector>

// =================================================================================================
//     ksw2 aligner
// =================================================================================================

// Two alignment variants:
//
//   align_ksw2_score — single forward EXTZ_ONLY|SCORE_ONLY pass.
//     Gives score + end position; start is not computed (has_start = false).
//
//   align_ksw2_cigar — reverse-first fitting alignment with CIGAR:
//     1. Reverse pass (SCORE_ONLY) on (rev_query, rev_target) → start_t = tlen-1-mqe_t
//     2. Forward pass (with CIGAR) on (query, target[start_t..]) → CIGAR + end_t
//     The forward pass is anchored at the known start, so the CIGAR is in the correct
//     forward orientation without reversal.
//
// Scoring: minimap2 short-read defaults — match=+2, mismatch=−4, gap-open=4, gap-extend=2.
// N scores 0 against any base (treated as wildcard).
// No hot/cold split — ksw2 has no precomputation step.

namespace {

inline uint8_t ksw2_encode( char c )
{
    switch( c ) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;
    }
}

// 5-symbol alphabet (ACGTN); rows/cols ordered A C G T N
constexpr int8_t KSW2_SCORE_MAT[25] = {
     2, -4, -4, -4,  0,
    -4,  2, -4, -4,  0,
    -4, -4,  2, -4,  0,
    -4, -4, -4,  2,  0,
     0,  0,  0,  0,  0
};

} // namespace

// Single forward pass: score + end position only (has_start = false).
inline AlignResult align_ksw2_score( std::string const& query, std::string const& target )
{
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    std::vector<uint8_t> q_enc( qlen ), t_enc( tlen );
    for( int i = 0; i < qlen; ++i ) q_enc[i] = ksw2_encode( query[i]  );
    for( int i = 0; i < tlen; ++i ) t_enc[i] = ksw2_encode( target[i] );

    int const flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY;

    ksw_extz_t ez = {};
#ifdef KSW2_HAS_SSE41
    ksw_extz2_sse(
        nullptr, qlen, q_enc.data(), tlen, t_enc.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, 0, flag, &ez
    );
#else
    ksw_extz(
        nullptr, qlen, q_enc.data(), tlen, t_enc.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, flag, &ez
    );
#endif

    if( ez.mqe == KSW_NEG_INF || ez.mqe_t < 0 ) {
        return AlignResult{ 0, 0, 0, true };
    }

    return AlignResult{
        0,              // start: not computed
        ez.mqe_t + 1,  // 0-based exclusive end in target
        ez.mqe,
        false,
        false           // has_start = false
    };
}

// Reverse-first two-pass alignment with CIGAR: score + start + end.
// Pass 1 (SCORE_ONLY): rev_query vs rev_target → start_t = tlen - 1 - mqe_t
// Pass 2 (with CIGAR): query vs target[start_t..] → end_t + CIGAR in forward orientation
inline AlignResult align_ksw2_cigar( std::string const& query, std::string const& target )
{
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    std::vector<uint8_t> q_enc( qlen ), t_enc( tlen );
    for( int i = 0; i < qlen; ++i ) q_enc[i] = ksw2_encode( query[i]  );
    for( int i = 0; i < tlen; ++i ) t_enc[i] = ksw2_encode( target[i] );

    // --- Reverse pass (SCORE_ONLY): find start position ---
    std::vector<uint8_t> rev_q( qlen ), rev_t( tlen );
    for( int i = 0; i < qlen; ++i ) rev_q[i] = q_enc[ qlen - 1 - i ];
    for( int i = 0; i < tlen; ++i ) rev_t[i] = t_enc[ tlen - 1 - i ];

    ksw_extz_t ez1 = {};
    int const flag_score = KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY;
#ifdef KSW2_HAS_SSE41
    ksw_extz2_sse(
        nullptr, qlen, rev_q.data(), tlen, rev_t.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, 0, flag_score, &ez1
    );
#else
    ksw_extz(
        nullptr, qlen, rev_q.data(), tlen, rev_t.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, flag_score, &ez1
    );
#endif

    if( ez1.mqe == KSW_NEG_INF || ez1.mqe_t < 0 ) {
        return AlignResult{ 0, 0, 0, true };
    }

    int const start_t   = tlen - 1 - ez1.mqe_t;  // 0-based inclusive start in target
    int const fwd_tlen  = tlen - start_t;

    // --- Forward pass (with CIGAR): anchored at start_t → end position + CIGAR ---
    ksw_extz_t ez2 = {};
#ifdef KSW2_HAS_SSE41
    ksw_extz2_sse(
        nullptr, qlen, q_enc.data(), fwd_tlen, t_enc.data() + start_t,
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, 0, KSW_EZ_EXTZ_ONLY, &ez2
    );
#else
    ksw_extz(
        nullptr, qlen, q_enc.data(), fwd_tlen, t_enc.data() + start_t,
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, KSW_EZ_EXTZ_ONLY, &ez2
    );
#endif

    if( ez2.mqe == KSW_NEG_INF || ez2.mqe_t < 0 ) {
        return AlignResult{ 0, 0, 0, true };
    }

    // CIGAR is in ez2.cigar (allocated by ksw2); free it — we don't store it in AlignResult.
    free( ez2.cigar );

    return AlignResult{
        start_t,                    // 0-based inclusive start in target
        start_t + ez2.mqe_t + 1,   // 0-based exclusive end in target
        ez2.mqe,
        false
    };
}
