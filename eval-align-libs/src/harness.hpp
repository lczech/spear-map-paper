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

#include "shatter.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>

// =================================================================================================
//     Filter
// =================================================================================================

// Returns true if the read is suitable for benchmarking.
// Rejects reads containing any ambiguous base (N), so that alignment errors are attributable
// to the aligner rather than the input quality.
inline bool passes_filter( ShatteredRead const& read )
{
    for( char const c : read.forward ) {
        if( c == 'N' ) {
            return false;
        }
    }
    return true;
}

// =================================================================================================
//     BenchmarkStats
// =================================================================================================

struct BenchmarkStats
{
    size_t total_reads    = 0;
    size_t filtered_reads = 0;
    size_t passing_reads  = 0;

    std::map<size_t, size_t> length_hist;   // read length (bp) → count
    std::map<size_t, size_t> window_hist;   // window size (bp) → count
};

// =================================================================================================
//     collect_stats
// =================================================================================================

// Streams through a shatter generator, applying the filter and collecting statistics.
// Templated on Range so it works with Generator<ShatteredRead> and any other compatible range.
template<typename Range>
BenchmarkStats collect_stats( Range&& reads )
{
    BenchmarkStats stats;

    for( auto const& read : reads ) {
        ++stats.total_reads;

        if( !passes_filter( read )) {
            ++stats.filtered_reads;
            continue;
        }
        ++stats.passing_reads;
        ++stats.length_hist[ read.origin_end - read.origin_start ];
        ++stats.window_hist[ read.window_end - read.window_start ];
    }

    return stats;
}

// =================================================================================================
//     print_stats
// =================================================================================================

inline void print_stats( BenchmarkStats const& stats )
{
    std::cout << "Total reads generated : " << stats.total_reads    << "\n";
    std::cout << "Filtered (contains N) : " << stats.filtered_reads << "\n";
    std::cout << "Passing reads         : " << stats.passing_reads  << "\n";

    // Read length histogram with ASCII bar scaled to 50 chars.
    // Deactivated for now to keep output succinct.
    // size_t const max_count = stats.length_hist.empty() ? 1 :
    //     std::max_element(
    //         stats.length_hist.begin(), stats.length_hist.end(),
    //         []( auto const& a, auto const& b ){ return a.second < b.second; }
    //     )->second;

    // std::cout << "\nRead length distribution:\n";
    // for( auto const& [len, count] : stats.length_hist ) {
    //     size_t const bar = ( count * 50 ) / max_count;
    //     std::cout
    //         << "  " << std::setw(4) << len << " bp | "
    //         << std::string( bar, '#' )
    //         << " " << count << "\n";
    // }

    // Window size histogram — typically just a few distinct values.
    std::cout << "\nReference window size distribution:\n";
    for( auto const& [size, count] : stats.window_hist ) {
        std::cout << "  " << std::setw(4) << size << " bp : " << count << "\n";
    }
}
