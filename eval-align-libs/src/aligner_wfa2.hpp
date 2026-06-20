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
#include "wfa.h"
#include "wavefront_align.h"
}

#include <string>

// =================================================================================================
//     WFA2 aligner
// =================================================================================================

// Four alignment variants (no hot/cold split — the aligner object is created once at startup
// and reused across all reads, modelling production use in Spear):
//
//   wfa2-score-exact    — compute_score, no heuristic; score + end only (has_start = false)
//   wfa2-score-heuristic— compute_score, wfadaptive(10,50,1); score + end only (has_start = false)
//   wfa2-cigar-exact    — compute_alignment, no heuristic; score + start + end
//   wfa2-cigar-heuristic— compute_alignment, wfadaptive(10,50,1); score + start + end
//
// Fitting (semi-global) alignment: pattern (query) fully consumed, text (reference) has free ends.
// text_begin_free and text_end_free are set to tlen per call via
// wavefront_aligner_set_alignment_free_ends(), since target length varies per read.
//
// Penalties: gap_affine, match=0, mismatch=4, gap_opening=4, gap_extension=2.
// Score in AlignResult is cigar->score: ≤0, with 0 = perfect alignment, more negative = worse.
// This is the opposite sign convention from a penalty (higher is better, same sign as parasail).
//
// Failure: any status other than WF_STATUS_ALG_COMPLETED sets failed=true.
// For heuristic variants, WF_STATUS_ALG_PARTIAL (heuristic terminated early) is also a failure.

// Factory: create a configured WFA2 aligner. Caller owns; free with wavefront_aligner_delete().
inline wavefront_aligner_t* make_wfa2_aligner( alignment_scope_t scope, bool heuristic )
{
    wavefront_aligner_attr_t attrs  = wavefront_aligner_attr_default;
    attrs.distance_metric           = gap_affine;
    attrs.alignment_scope           = scope;
    attrs.alignment_form.span       = alignment_endsfree;
    // pattern (query) consumed fully; text (reference) free ends set per-call
    attrs.alignment_form.pattern_begin_free = 0;
    attrs.alignment_form.pattern_end_free   = 0;
    attrs.alignment_form.text_begin_free    = 0;
    attrs.alignment_form.text_end_free      = 0;
    attrs.affine_penalties.match         = 0;
    attrs.affine_penalties.mismatch      = 4;
    attrs.affine_penalties.gap_opening   = 4;
    attrs.affine_penalties.gap_extension = 2;
    attrs.system.verbose = 0;

    wavefront_aligner_t* wf = wavefront_aligner_new( &attrs );

    if( heuristic ) {
        wavefront_aligner_set_heuristic_wfadaptive( wf, 10, 50, 1 );
    }

    return wf;
}

// Score-only: single pass, score + end position (has_start = false).
inline AlignResult align_wfa2_score(
    wavefront_aligner_t* wf,
    std::string const& query,
    std::string const& target
) {
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    wavefront_aligner_set_alignment_free_ends( wf, 0, 0, tlen, tlen );
    int const status = wavefront_align( wf, query.c_str(), qlen, target.c_str(), tlen );

    if( status != WF_STATUS_ALG_COMPLETED ) {
        return AlignResult{ 0, 0, 0, true };
    }

    // In compute_score mode, wavefront_unialign_terminate sets cig->end_h = text_length
    // (discarding the actual fitting endpoint). The real end is in alignment_end_pos.offset,
    // which wavefront_termination_endsfree sets before the terminate call.
    return AlignResult{
        0,                                  // start: not computed
        wf->alignment_end_pos.offset,       // 0-based exclusive end in target
        wf->cigar->score,                   // WFA2 penalty (≥0, lower = better)
        false,
        false                               // has_start = false
    };
}

// Full alignment: score + start + end, start derived by counting reference-consuming CIGAR ops.
inline AlignResult align_wfa2_cigar(
    wavefront_aligner_t* wf,
    std::string const& query,
    std::string const& target
) {
    int const qlen = static_cast<int>( query.size() );
    int const tlen = static_cast<int>( target.size() );

    wavefront_aligner_set_alignment_free_ends( wf, 0, 0, tlen, tlen );
    int const status = wavefront_align( wf, query.c_str(), qlen, target.c_str(), tlen );

    if( status != WF_STATUS_ALG_COMPLETED ) {
        return AlignResult{ 0, 0, 0, true };
    }

    cigar_t* const cig = wf->cigar;

    // Count reference-consuming ops (M, X, D) to derive start position.
    // With has_misms=false (default), both matches and mismatches are encoded as 'M'.
    int ref_span = 0;
    for( int i = cig->begin_offset; i < cig->end_offset; ++i ) {
        char const op = cig->operations[i];
        if( op == 'M' || op == 'X' || op == 'D' ) ++ref_span;
    }

    return AlignResult{
        cig->end_h - ref_span,  // 0-based inclusive start in target
        cig->end_h,             // 0-based exclusive end in target
        cig->score,
        false
    };
}
