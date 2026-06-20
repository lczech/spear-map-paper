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

#include "aligner_edlib.hpp"
#include "aligner_ksw2.hpp"
#include "aligner_parasail.hpp"
#include "aligner_wfa2.hpp"
#include "stats.hpp"
#include "mutate.hpp"
#include "shatter.hpp"

#include "genesis/util/core/logging.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

// =================================================================================================
//      Main Program
// =================================================================================================

int main( int argc, char** argv )
{
    // -------------------------------------------------------------------------
    //     Init
    // -------------------------------------------------------------------------

    genesis::util::core::Logging::log_to_stdout();
    genesis::util::core::Logging::details.level = false;
    genesis::util::core::Logging::details.time  = true;
    genesis::util::core::Options::get().init_global_thread_pool( 1 );

    // -------------------------------------------------------------------------
    //     Arguments
    // -------------------------------------------------------------------------

    std::string fasta_path;
    size_t max_reads = 0;         // 0 = unlimited
    size_t max_chromosomes = 0;   // 0 = unlimited

    for( int a = 1; a < argc; ++a ) {
        std::string const arg = argv[a];
        if(( arg == "--max-reads" || arg == "--max-chromosomes" ) && a + 1 < argc ) {
            size_t const val = std::stoull( argv[++a] );
            if( arg == "--max-reads" )       max_reads       = val;
            if( arg == "--max-chromosomes" ) max_chromosomes = val;
        } else if( fasta_path.empty() && arg[0] != '-' ) {
            fasta_path = arg;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            std::cerr << "Usage: " << argv[0]
                      << " <reference.fasta>"
                      << " [--max-reads N]"
                      << " [--max-chromosomes N]\n";
            return 1;
        }
    }
    if( fasta_path.empty() ) {
        std::cerr << "Usage: " << argv[0]
                  << " <reference.fasta>"
                  << " [--max-reads N]"
                  << " [--max-chromosomes N]\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    //     Mutation parameter grid
    // -------------------------------------------------------------------------

    // Each entry defines one combination of divergence + damage to benchmark.
    // sub_rate and indel_rate model evolutionary divergence from the reference.
    // damage_rate models aDNA deamination (C→T at 5', G→A at 3').
    std::vector<MutateParams> const grid = {
        // sub_rate indel_rate indel_mean_len damage_rate decay_lambda
        {  0.00,    0.00,      1.5,           0.00,       0.3  }, // clean baseline
        {  0.05,    0.00,      1.5,           0.00,       0.3  }, // 5% divergence only
        {  0.10,    0.00,      1.5,           0.00,       0.3  }, // 10% divergence only
        {  0.05,    0.02,      1.5,           0.00,       0.3  }, // div + indels
        {  0.05,    0.00,      1.5,           0.10,       0.3  }, // div + damage
        {  0.05,    0.02,      1.5,           0.10,       0.3  }, // div + indels + damage
    };

    // -------------------------------------------------------------------------
    //     Shatter and collect stats
    // -------------------------------------------------------------------------

    LOG_MSG << "Alignment library benchmarking";

    // Parasail matrices initialised once; custom is mutable (freed below), dnafull is built-in.
    parasail_matrix_t*       mat_custom  = parasail_make_custom_matrix();
    parasail_matrix_t const* mat_dnafull = parasail_matrix_lookup( "dnafull" );

    // WFA2 aligners created once and reused across all reads (models production use in Spear).
    wavefront_aligner_t* wf_score_exact     = make_wfa2_aligner( compute_score,     false );
    wavefront_aligner_t* wf_score_heuristic = make_wfa2_aligner( compute_score,     true  );
    wavefront_aligner_t* wf_cigar_exact     = make_wfa2_aligner( compute_alignment, false );
    wavefront_aligner_t* wf_cigar_heuristic = make_wfa2_aligner( compute_alignment, true  );

    std::mt19937_64 rng( 42 );
    std::vector<std::vector<size_t>> length_hists( grid.size() );

    // One map per grid entry: aligner-name → stats.
    // Map iteration is alphabetical, so names are chosen to give a sensible print order.
    // Variants we benchmark:
    //   edlib        — edit-distance HW alignment; baseline, no affine gaps
    //   ksw2         — affine-gap HW alignment via reverse trick (two extension passes)
    //
    // parasail variants test three orthogonal axes:
    //   Matrix  : custom (match=+2 / mismatch=-4 / N=0, same as ksw2 for direct comparison)
    //             vs dnafull (IUPAC-aware built-in matrix, for real-data N handling)
    //   Approach: score = single forward pass, score + end only (has_start = false)
    //             cigar = single pass with traceback, gives CIGAR and exact start/end directly
    //   Timing  : cold = profile creation + alignment (first time seeing this query)
    //             hot  = alignment only (profile already built, models repeated use)
    std::vector<std::map<std::string, BenchmarkStats>> stats( grid.size() );
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( auto const& name : {
            "edlib",
            "ksw2-cigar",
            "ksw2-score",
            "parasail-cigar-custom-cold",
            "parasail-cigar-custom-hot",
            "parasail-cigar-dnafull-cold",
            "parasail-cigar-dnafull-hot",
            "parasail-score-custom-cold",
            "parasail-score-custom-hot",
            "parasail-score-dnafull-cold",
            "parasail-score-dnafull-hot",
            "wfa2-cigar-exact",
            "wfa2-cigar-heuristic",
            "wfa2-score-exact",
            "wfa2-score-heuristic",
        }) {
            stats[i].emplace( name, BenchmarkStats( name ) );
        }
    }

    size_t total_reads      = 0;
    size_t filtered_reads   = 0;
    size_t chromosomes_seen = 0;
    std::string last_chromosome;
    std::map<size_t, size_t> window_hist;

    for( auto const& read : shatter( fasta_path )) {
        // Limits for quick testing
        if( read.chromosome != last_chromosome ) {
            last_chromosome = read.chromosome;
            ++chromosomes_seen;
            if( max_chromosomes > 0 && chromosomes_seen > max_chromosomes ) {
                break;
            }
        }
        ++total_reads;
        if( max_reads > 0 && total_reads > max_reads ) {
            break;
        }

        // Filter out anything with Ns - enough for benchmarking.
        if( !passes_filter( read )) {
            ++filtered_reads;
            continue;
        }
        ++window_hist[ read.window_end - read.window_start ];

        // True read position within the reference window (same for all param sets)
        int32_t const true_start = static_cast<int32_t>( read.origin_start - read.window_start );
        int32_t const true_end   = static_cast<int32_t>( read.origin_end   - read.window_start );

        // Mutate and align for each parameter set in the grid of mutation params.
        for( size_t i = 0; i < grid.size(); ++i ) {
            auto const mutated = mutate( read, grid[i], rng );

            size_t const len = mutated.forward.size();
            if( len >= length_hists[i].size() ) length_hists[i].resize( len + 1, 0 );
            ++length_hists[i][len];

            // Helper lambda to time a single alignment call and accumulate its result.
            auto const time_align = [&]( auto align_fn, std::string const& key ) {
                auto const t0 = std::chrono::high_resolution_clock::now();
                AlignResult const result = align_fn();
                auto const t1 = std::chrono::high_resolution_clock::now();
                uint64_t const ns = static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>( t1 - t0 ).count()
                );
                accumulate_align_result( stats[i][key], result, ns, true_start, true_end );
            };

            // --- edlib (no hot/cold split — no precomputation) ---
            time_align( [&]{
                return align_edlib( mutated.forward, mutated.window );
            }, "edlib" );

            // --- ksw2 ---
            time_align( [&]{
                return align_ksw2_score( mutated.forward, mutated.window );
            }, "ksw2-score" );
            time_align( [&]{
                return align_ksw2_cigar( mutated.forward, mutated.window );
            }, "ksw2-cigar" );

            // --- parasail ---
            // score-custom: single forward pass, custom scoring matrix (match=+2/mismatch=-4/N=0)
            {
                // hot: profile pre-built, only alignment timed
                auto* pf = parasail_make_profile( mutated.forward, mat_custom );
                time_align( [&]{
                    return align_parasail_score( pf, mutated.window );
                }, "parasail-score-custom-hot" );
                parasail_profile_free( pf );

                // cold: profile creation + alignment timed together
                time_align( [&]{
                    auto* p = parasail_make_profile( mutated.forward, mat_custom );
                    auto r  = align_parasail_score( p, mutated.window );
                    parasail_profile_free( p );
                    return r;
                }, "parasail-score-custom-cold" );
            }

            // score-dnafull: single forward pass, IUPAC-aware DNAfull matrix
            {
                auto* pf = parasail_make_profile( mutated.forward, mat_dnafull );
                time_align( [&]{ return align_parasail_score( pf, mutated.window ); },
                         "parasail-score-dnafull-hot" );
                parasail_profile_free( pf );

                time_align( [&]{
                    auto* p = parasail_make_profile( mutated.forward, mat_dnafull );
                    auto r  = align_parasail_score( p, mutated.window );
                    parasail_profile_free( p );
                    return r;
                }, "parasail-score-dnafull-cold" );
            }

            // cigar-custom: single traceback pass, beg_ref from cigar struct, custom matrix
            {
                auto* pf = parasail_make_profile( mutated.forward, mat_custom );
                time_align( [&]{ return align_parasail_cigar( pf, mutated.forward, mutated.window, mat_custom ); },
                         "parasail-cigar-custom-hot" );
                parasail_profile_free( pf );

                time_align( [&]{
                    auto* p = parasail_make_profile( mutated.forward, mat_custom );
                    auto r  = align_parasail_cigar( p, mutated.forward, mutated.window, mat_custom );
                    parasail_profile_free( p );
                    return r;
                }, "parasail-cigar-custom-cold" );
            }

            // cigar-dnafull: single traceback pass, DNAfull matrix
            {
                auto* pf = parasail_make_profile( mutated.forward, mat_dnafull );
                time_align( [&]{ return align_parasail_cigar( pf, mutated.forward, mutated.window, mat_dnafull ); },
                         "parasail-cigar-dnafull-hot" );
                parasail_profile_free( pf );

                time_align( [&]{
                    auto* p = parasail_make_profile( mutated.forward, mat_dnafull );
                    auto r  = align_parasail_cigar( p, mutated.forward, mutated.window, mat_dnafull );
                    parasail_profile_free( p );
                    return r;
                }, "parasail-cigar-dnafull-cold" );
            }

            // --- wfa2 ---
            time_align( [&]{
                return align_wfa2_score( wf_score_exact, mutated.forward, mutated.window );
            }, "wfa2-score-exact" );
            time_align( [&]{
                return align_wfa2_score( wf_score_heuristic, mutated.forward, mutated.window );
            }, "wfa2-score-heuristic" );
            time_align( [&]{
                return align_wfa2_cigar( wf_cigar_exact, mutated.forward, mutated.window );
            }, "wfa2-cigar-exact" );
            time_align( [&]{
                return align_wfa2_cigar( wf_cigar_heuristic, mutated.forward, mutated.window );
            }, "wfa2-cigar-heuristic" );

        }
    }
    LOG_MSG << "Finished processing";

    // -------------------------------------------------------------------------
    //     Report
    // -------------------------------------------------------------------------

    print_global_stats( total_reads, filtered_reads, window_hist );

    for( size_t i = 0; i < grid.size(); ++i ) {
        print_param_header( grid[i], length_hists[i] );
        for( auto const& [name, s] : stats[i] ) {
            print_stats( s );
        }
    }

    write_all_csvs( BENCH_DIR, grid, length_hists, stats, window_hist );

    parasail_matrix_free( mat_custom );
    wavefront_aligner_delete( wf_score_exact     );
    wavefront_aligner_delete( wf_score_heuristic );
    wavefront_aligner_delete( wf_cigar_exact     );
    wavefront_aligner_delete( wf_cigar_heuristic );
    return 0;
}
