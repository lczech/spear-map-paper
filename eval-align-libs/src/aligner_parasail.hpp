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
#include "parasail.h"
}

#include <string>

// =================================================================================================
//     parasail aligner
// =================================================================================================

// Fitting (semi-global) alignment: parasail sg_dx mode — free gaps at the start and end of the
// reference, query must be fully consumed. This is the same mode as edlib's EDLIB_MODE_HW;
// "HW" is edlib-specific naming, the general term is fitting or infix alignment.
//
// Scoring: gap-open=4, gap-extend=2 (minimap2 short-read defaults, matching ksw2).
// SIMD width and fallback selected automatically at runtime via _sat (8-bit with 16-bit retry).
//
// Two alignment variants:
//   loc   — reverse trick (two sg_dx passes), no CIGAR, start derived from second pass
//   cigar — single sg_dx_trace pass; beg_ref/end_ref read from parasail_cigar_t directly
//
// Each function takes pre-built profile(s).  Cold/hot timing split is handled in main:
//   cold = profile creation + alignment (timed together)
//   hot  = alignment only (profiles pre-built outside the timer)

// =================================================================================================
//     Matrix helpers
// =================================================================================================

// 5-symbol ACGTN matrix: match=+2, mismatch=−4, N scores 0 against any base.
// Caller owns the returned pointer; free with parasail_matrix_free().
inline parasail_matrix_t* parasail_make_custom_matrix()
{
    parasail_matrix_t* mat = parasail_matrix_create( "ACGTN", 2, -4 );
    for( int j = 0; j < 5; ++j ) {
        parasail_matrix_set_value( mat, 4, j, 0 );
        parasail_matrix_set_value( mat, j, 4, 0 );
    }
    return mat;
}

// Convenience wrapper: build a sat profile from a string and a matrix.
// Caller owns the returned pointer; free with parasail_profile_free().
inline parasail_profile_t* parasail_make_profile(
    std::string const& seq,
    parasail_matrix_t const* matrix
) {
    return parasail_profile_create_sat(
        seq.c_str(), static_cast<int>( seq.size() ), matrix
    );
}

// =================================================================================================
//     Location-only alignment (reverse trick, no CIGAR)
// =================================================================================================

// pf = profile of the forward query, pr = profile of the reversed query.
// Forward pass → end_ref; reverse pass on rev(target[0..end_ref]) → start_ref.
inline AlignResult align_parasail_loc(
    parasail_profile_t* pf,
    parasail_profile_t* pr,
    std::string const& target
) {
    int const tlen = static_cast<int>( target.size() );

    // Forward pass
    parasail_result_t* rf = parasail_sg_dx_striped_profile_sat(
        pf, target.c_str(), tlen, 4, 2
    );
    if( !rf ) return AlignResult{ 0, 0, 0, true };

    int const end_t = rf->end_ref;
    int const score = rf->score;
    parasail_result_free( rf );

    if( end_t < 0 ) return AlignResult{ 0, 0, 0, true };

    // Reverse pass: rev(query) vs reverse(target[0..end_t])
    int const rev_tlen = end_t + 1;
    std::string const rev_target(
        target.crbegin() + ( tlen - 1 - end_t ),
        target.crend()
    );

    parasail_result_t* rr = parasail_sg_dx_striped_profile_sat(
        pr, rev_target.c_str(), rev_tlen, 4, 2
    );
    int const start_t = ( rr && rr->end_ref >= 0 ) ? ( end_t - rr->end_ref ) : 0;
    if( rr ) parasail_result_free( rr );

    return AlignResult{
        start_t,      // 0-based inclusive start in target
        end_t + 1,    // 0-based exclusive end   in target (edlib convention)
        score,
        false
    };
}

// =================================================================================================
//     CIGAR alignment (single pass with traceback)
// =================================================================================================

// pf = forward query profile. query and target strings are needed for parasail_result_get_cigar.
// beg_ref is read directly from the cigar struct — no second pass required.
inline AlignResult align_parasail_cigar(
    parasail_profile_t* pf,
    std::string const& query,
    std::string const& target,
    parasail_matrix_t const* matrix
) {
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    parasail_result_t* r = parasail_sg_dx_trace_striped_profile_sat(
        pf, target.c_str(), tlen, 4, 2
    );
    if( !r ) return AlignResult{ 0, 0, 0, true };

    int const end_t = r->end_ref;
    int const score = r->score;

    parasail_cigar_t* cigar = parasail_result_get_cigar(
        r, query.c_str(), qlen, target.c_str(), tlen, matrix
    );
    int const start_t = cigar ? cigar->beg_ref : 0;
    if( cigar ) parasail_cigar_free( cigar );
    parasail_result_free( r );

    return AlignResult{
        start_t,
        end_t + 1,
        score,
        false
    };
}
