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
#include "stats.hpp"
#include "mutate.hpp"
#include "shatter.hpp"

#include "genesis/util/core/logging.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
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

    if( argc != 2 ) {
        std::cerr << "Usage: " << argv[0] << " <reference.fasta>\n";
        return 1;
    }
    std::string const fasta_path = argv[1];

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

    std::mt19937_64 rng( 42 );
    std::vector<BenchmarkStats> stats;
    stats.reserve( grid.size() );
    for( size_t i = 0; i < grid.size(); ++i ) {
        stats.emplace_back( "edlib" );
    }

    size_t total_reads    = 0;
    size_t filtered_reads = 0;
    std::map<size_t, size_t> window_hist;

    for( auto const& read : shatter( fasta_path )) {
        ++total_reads;
        if( !passes_filter( read )) {
            ++filtered_reads;
            continue;
        }
        ++window_hist[ read.window_end - read.window_start ];

        // True read position within the reference window (same for all param sets)
        int32_t const true_start = static_cast<int32_t>( read.origin_start - read.window_start );
        int32_t const true_end   = static_cast<int32_t>( read.origin_end   - read.window_start );

        for( size_t i = 0; i < grid.size(); ++i ) {
            auto const mutated = mutate( read, grid[i], rng );
            ++stats[i].passing_reads;
            ++stats[i].length_hist[ mutated.forward.size() ];

            // --- edlib alignment (no hot/cold split — edlib has no precomputation) ---
            auto const t0 = std::chrono::high_resolution_clock::now();
            AlignResult const result = align_edlib( mutated.forward, mutated.window );
            auto const t1 = std::chrono::high_resolution_clock::now();
            uint64_t const ns = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>( t1 - t0 ).count()
            );
            accumulate_align_result( stats[i], result, ns, true_start, true_end );
        }
    }

    // -------------------------------------------------------------------------
    //     Report
    // -------------------------------------------------------------------------

    print_global_stats( total_reads, filtered_reads, window_hist );

    for( size_t i = 0; i < grid.size(); ++i ) {
        print_stats( grid[i], stats[i] );
    }

    return 0;
}
