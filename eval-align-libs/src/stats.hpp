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

#include "mutate.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>

// =================================================================================================
//     BenchmarkStats
// =================================================================================================

// Per-mutation-parameter-set statistics.
// Window size histogram and total/filtered counts are global and tracked separately in main.
struct BenchmarkStats
{
    size_t passing_reads = 0;
    std::map<size_t, size_t> length_hist;   // mutated read length (bp) → count
};

// =================================================================================================
//     print_stats
// =================================================================================================

inline void print_global_stats(
    size_t total_reads,
    size_t filtered_reads,
    std::map<size_t, size_t> const& window_hist
) {
    std::cout << "Total reads generated : " << total_reads    << "\n";
    std::cout << "Filtered (contains N) : " << filtered_reads << "\n";
    std::cout << "Passing reads         : " << ( total_reads - filtered_reads ) << "\n";

    std::cout << "\nReference window size distribution:\n";
    for( auto const& [size, count] : window_hist ) {
        std::cout << "  " << std::setw(4) << size << " bp : " << count << "\n";
    }
}

inline void print_stats( MutateParams const& params, BenchmarkStats const& stats )
{
    std::cout << std::fixed << std::setprecision( 3 );
    std::cout
        << "\n[ sub=" << params.sub_rate
        << "  indel=" << params.indel_rate
        << "  mean_len=" << params.indel_mean_len
        << "  dmg=" << params.damage_rate
        << "  lambda=" << params.decay_lambda
        << " ]  passing=" << stats.passing_reads << "\n";

    if( stats.length_hist.empty() ) {
        return;
    }

    size_t const max_count = std::max_element(
        stats.length_hist.begin(), stats.length_hist.end(),
        []( auto const& a, auto const& b ){ return a.second < b.second; }
    )->second;

    for( auto const& [len, count] : stats.length_hist ) {
        size_t const bar = ( count * 40 ) / max_count;
        std::cout
            << "  " << std::setw(4) << len << " bp | "
            << std::string( bar, '#' ) << " " << count << "\n";
    }
}
