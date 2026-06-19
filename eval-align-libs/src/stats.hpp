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
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

// =================================================================================================
//     AlignResult
// =================================================================================================

// Uniform return type for all aligner wrappers.
// start and end are 0-based positions within the reference window.
// score is whatever the library natively returns (edit distance, alignment score, etc.).
struct AlignResult
{
    int32_t start  = 0;
    int32_t end    = 0;
    int32_t score  = 0;
    bool    failed = false;
};

// =================================================================================================
//     BenchmarkStats
// =================================================================================================

// Per-aligner-call statistics (one instance per library × hot/cold variant).
// Name identifies the aligner and call type (e.g. "edlib", "ksw2-cold", "ksw2-hot").
// Window size histogram and total/filtered counts are global and tracked separately in main.
struct BenchmarkStats
{
    explicit BenchmarkStats( std::string n = {} ) : name( std::move(n) ) {}

    std::string name;

    size_t successful_alignments = 0;
    size_t failed_alignments     = 0;

    uint64_t total_ns = 0;
    uint64_t min_ns   = std::numeric_limits<uint64_t>::max();
    uint64_t max_ns   = 0;

    std::vector<uint64_t> timing_hist;            // log10 buckets: 0=<1µs, 1=1-10µs, …

    std::map<int32_t, size_t> start_offset_hist; // aligner start - true start → count
    std::map<int32_t, size_t> end_offset_hist;   // aligner end   - true end   → count
    std::map<int32_t, size_t> score_hist;        // native library score/edit distance → count
};

// =================================================================================================
//     Timing histogram helpers
// =================================================================================================

// Bucket 0: < 1 µs (ns < 1000)
// Bucket n (n >= 1): [10^(n+2), 10^(n+3)) ns, i.e. 1-10 µs, 10-100 µs, 100 µs-1 ms, …
inline size_t timing_bucket( uint64_t ns )
{
    if( ns < 1000 ) return 0;
    return static_cast<size_t>( std::log10( static_cast<double>( ns ))) - 2;
}

inline std::string timing_bucket_label( size_t bucket )
{
    if( bucket == 0 ) return "<1 us";

    uint64_t lo = 1;
    for( size_t i = 0; i < bucket + 2; ++i ) lo *= 10;
    uint64_t const hi = lo * 10;

    auto fmt = []( uint64_t ns ) -> std::string {
        if( ns < 1'000'000ull )     return std::to_string( ns / 1'000ull )         + " us";
        if( ns < 1'000'000'000ull ) return std::to_string( ns / 1'000'000ull )     + " ms";
        return                             std::to_string( ns / 1'000'000'000ull ) + " s";
    };
    return fmt(lo) + "-" + fmt(hi);
}

// =================================================================================================
//     accumulate_align_result
// =================================================================================================

inline void accumulate_align_result(
    BenchmarkStats& stats,
    AlignResult const& result,
    uint64_t ns,
    int32_t true_start_in_window,
    int32_t true_end_in_window
) {
    stats.total_ns += ns;
    stats.min_ns    = std::min( stats.min_ns, ns );
    stats.max_ns    = std::max( stats.max_ns, ns );

    size_t const b = timing_bucket( ns );
    if( b >= stats.timing_hist.size() ) {
        stats.timing_hist.resize( b + 1, 0 );
    }
    ++stats.timing_hist[b];

    if( result.failed ) {
        ++stats.failed_alignments;
        return;
    }

    ++stats.successful_alignments;
    ++stats.start_offset_hist[ result.start - true_start_in_window ];
    ++stats.end_offset_hist[   result.end   - true_end_in_window   ];
    ++stats.score_hist[ result.score ];
}

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

inline void print_hist(
    std::string const& label,
    std::map<int32_t, size_t> const& hist,
    std::string const& unit = ""
) {
    if( hist.empty() ) return;
    size_t const max_count = std::max_element(
        hist.begin(), hist.end(),
        []( auto const& a, auto const& b ){ return a.second < b.second; }
    )->second;
    std::cout << "  " << label << ":\n";
    for( auto const& [val, count] : hist ) {
        size_t const bar = ( count * 40 ) / max_count;
        std::cout
            << "    " << std::setw(5) << val << unit << " | "
            << std::string( bar, '#' ) << " " << count << "\n";
    }
}

inline void print_timing_hist( std::vector<uint64_t> const& hist )
{
    if( hist.empty() ) return;
    uint64_t const max_count = *std::max_element( hist.begin(), hist.end() );
    if( max_count == 0 ) return;

    std::cout << "  timing:\n";
    for( size_t b = 0; b < hist.size(); ++b ) {
        size_t const bar = static_cast<size_t>( ( hist[b] * 40 ) / max_count );
        std::cout
            << "    " << std::setw(14) << timing_bucket_label(b) << " | "
            << std::string( bar, '#' ) << " " << hist[b] << "\n";
    }
}

inline void print_length_hist( std::vector<size_t> const& hist )
{
    if( hist.empty() ) return;
    size_t const max_count = *std::max_element( hist.begin(), hist.end() );
    if( max_count == 0 ) return;

    std::cout << "  read length (bp):\n";
    for( size_t len = 0; len < hist.size(); ++len ) {
        if( hist[len] == 0 ) continue;
        size_t const bar = ( hist[len] * 40 ) / max_count;
        std::cout
            << "    " << std::setw(5) << len << " bp | "
            << std::string( bar, '#' ) << " " << hist[len] << "\n";
    }
}

inline void print_param_header( MutateParams const& params, std::vector<size_t> const& length_hist )
{
    std::cout << "\n###########################################################################\n";
    std::cout << std::fixed << std::setprecision( 3 );
    std::cout
        << "sub=" << params.sub_rate
        << "  indel=" << params.indel_rate
        << "  mean_len=" << params.indel_mean_len
        << "  dmg=" << params.damage_rate
        << "  lambda=" << params.decay_lambda
        << "\n";
    print_length_hist( length_hist );
}

inline void print_stats( BenchmarkStats const& stats )
{
    std::cout << "\n---------------------------------------------------------------------------\n";
    std::cout << "[ " << stats.name << " ]\n";
    size_t const total = stats.successful_alignments + stats.failed_alignments;
    std::cout << "  processed=" << total;
    std::cout << "  successful=" << stats.successful_alignments;
    std::cout << "  failed=" << stats.failed_alignments;
    std::cout << "\n";

    if( stats.total_ns > 0 && total > 0 ) {
        std::cout << std::fixed << std::setprecision( 1 );
        std::cout
            << "  timing: mean=" << ( stats.total_ns / total )
            << " ns  min=" << stats.min_ns
            << " ns  max=" << stats.max_ns << " ns\n";
        print_timing_hist( stats.timing_hist );
    }

    print_hist( "start offset", stats.start_offset_hist, " bp" );
    print_hist( "end offset",   stats.end_offset_hist,   " bp" );
    print_hist( "score",        stats.score_hist );
}
