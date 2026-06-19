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

// Fitting (semi-global) alignment via the reverse trick:
//   1. Forward pass on (query, target)          → mqe_t = inclusive end in target
//   2. Reverse pass on (rev_query, rev_target[0..mqe_t]) → mqe_t gives start position
// Both passes use EXTZ_ONLY | SCORE_ONLY (no CIGAR, no zigzag DP).
// Scoring: minimap2 short-read defaults — match=+2, mismatch=−4, gap-open=4, gap-extend=2.
// N scores 0 against any base (treated as wildcard).

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

inline AlignResult align_ksw2( std::string const& query, std::string const& target )
{
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    std::vector<uint8_t> q_enc( qlen ), t_enc( tlen );
    for( int i = 0; i < qlen; ++i ) q_enc[i] = ksw2_encode( query[i]  );
    for( int i = 0; i < tlen; ++i ) t_enc[i] = ksw2_encode( target[i] );

    int const flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY;

    // --- Forward pass: query vs full target → end position ---
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

    int const end_t = ez.mqe_t;   // 0-based inclusive end in target
    int const score = ez.mqe;

    // --- Reverse pass: rev(query) vs rev(target[0..end_t]) → start position ---
    int const rev_tlen = end_t + 1;
    std::vector<uint8_t> rev_q( qlen ), rev_t( rev_tlen );
    for( int i = 0; i < qlen;     ++i ) rev_q[i] = q_enc[ qlen  - 1 - i ];
    for( int i = 0; i < rev_tlen; ++i ) rev_t[i] = t_enc[ end_t     - i ];

    ksw_extz_t ez2 = {};
#ifdef KSW2_HAS_SSE41
    ksw_extz2_sse(
        nullptr, qlen, rev_q.data(), rev_tlen, rev_t.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, 0, flag, &ez2
    );
#else
    ksw_extz(
        nullptr, qlen, rev_q.data(), rev_tlen, rev_t.data(),
        5, KSW2_SCORE_MAT, 4, 2, -1, -1, flag, &ez2
    );
#endif

    int const start_t = ( ez2.mqe_t >= 0 ) ? ( end_t - ez2.mqe_t ) : 0;

    return AlignResult{
        start_t,      // 0-based inclusive start in target
        end_t + 1,    // 0-based exclusive end  in target (edlib convention)
        score,
        false
    };
}
