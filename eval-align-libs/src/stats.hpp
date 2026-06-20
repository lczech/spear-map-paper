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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
    int32_t start     = 0;
    int32_t end       = 0;
    int32_t score     = 0;
    bool    failed    = false;
    bool    has_start = true;   // false for score-only variants (start is not computed)
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

    std::vector<uint64_t> timing_hist;            // 10 per-decade log10 buckets (see timing_bucket)

    std::map<int32_t, size_t> start_offset_hist; // aligner start - true start → count
    std::map<int32_t, size_t> end_offset_hist;   // aligner end   - true end   → count
    std::map<int32_t, size_t> score_hist;        // native library score/edit distance → count
};

// =================================================================================================
//     Timing histogram helpers
// =================================================================================================

// 10 uniform sub-buckets per decade on a log10 scale, starting at 1 ns.
// Bucket 0 is reserved for ns=0 (should not occur in practice).
// Bucket b covers [10^((b-1)/10), 10^(b/10)) ns  for b >= 1.
inline size_t timing_bucket( uint64_t ns )
{
    if( ns == 0 ) return 0;
    return static_cast<size_t>( std::floor( 10.0 * std::log10( static_cast<double>(ns) ) ) ) + 1;
}

// Returns [low_ns, high_ns) for the given bucket index (approximate integer bounds).
inline std::pair<uint64_t, uint64_t> timing_bucket_bounds( size_t bucket )
{
    if( bucket == 0 ) return { 0, 1 };
    size_t const b  = bucket - 1;
    uint64_t const lo = static_cast<uint64_t>( std::pow( 10.0, b / 10.0 ) );
    uint64_t const hi = static_cast<uint64_t>( std::ceil( std::pow( 10.0, (b + 1) / 10.0 ) ) );
    return { lo, std::max( hi, lo + 1 ) };
}

inline std::string timing_bucket_label( size_t bucket )
{
    if( bucket == 0 ) return "<1 us";
    auto [lo, hi] = timing_bucket_bounds( bucket );
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
    if( result.has_start ) {
        ++stats.start_offset_hist[ result.start - true_start_in_window ];
    }
    ++stats.end_offset_hist[ result.end - true_end_in_window ];
    ++stats.score_hist[ result.score ];
}

// =================================================================================================
//     stdout printing  (histograms commented out; only scalars remain)
// =================================================================================================

inline void print_global_stats(
    size_t total_reads,
    size_t filtered_reads,
    std::map<size_t, size_t> const& /* window_hist */
) {
    std::cout << "Total reads generated : " << total_reads    << "\n";
    std::cout << "Filtered (contains N) : " << filtered_reads << "\n";
    std::cout << "Passing reads         : " << ( total_reads - filtered_reads ) << "\n";

    // Window size distribution — written to windows.csv instead
    // std::cout << "\nReference window size distribution:\n";
    // for( auto const& [size, count] : window_hist ) {
    //     std::cout << "  " << std::setw(4) << size << " bp : " << count << "\n";
    // }
}

// Histogram bar printers kept for reference but no longer called from print_stats.
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

inline void print_param_header( MutateParams const& params, std::vector<size_t> const& /* length_hist */ )
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
    // Length histogram — written to lengths.csv instead
    // print_length_hist( length_hist );
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
        // Timing and offset histograms — written to CSV files instead
        // print_timing_hist( stats.timing_hist );
    }

    // print_hist( "start offset", stats.start_offset_hist, " bp" );
    // print_hist( "end offset",   stats.end_offset_hist,   " bp" );
    // print_hist( "score",        stats.score_hist );
}

// =================================================================================================
//     CSV output
// =================================================================================================

// Helper: write the grid-parameter prefix columns to an open CSV row.
inline void csv_grid_prefix( std::ofstream& f, size_t grid_idx, MutateParams const& p )
{
    f << grid_idx << ","
      << std::fixed << std::setprecision( 3 )
      << p.sub_rate       << ","
      << p.indel_rate     << ","
      << p.indel_mean_len << ","
      << p.damage_rate    << ","
      << p.decay_lambda   << ",";
}

inline void write_lengths_csv(
    std::string const& path,
    std::vector<MutateParams> const& grid,
    std::vector<std::vector<size_t>> const& length_hists
) {
    std::ofstream f( path );
    f << "grid_idx,sub_rate,indel_rate,indel_mean_len,damage_rate,decay_lambda,length_bp,count\n";
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( size_t len = 0; len < length_hists[i].size(); ++len ) {
            if( length_hists[i][len] == 0 ) continue;
            csv_grid_prefix( f, i, grid[i] );
            f << len << "," << length_hists[i][len] << "\n";
        }
    }
}

inline void write_timing_csv(
    std::string const& path,
    std::vector<MutateParams> const& grid,
    std::vector<std::map<std::string, BenchmarkStats>> const& stats
) {
    std::ofstream f( path );
    f << "grid_idx,sub_rate,indel_rate,indel_mean_len,damage_rate,decay_lambda,"
         "aligner,bucket_low_ns,bucket_high_ns,count\n";
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( auto const& [name, s] : stats[i] ) {
            for( size_t b = 0; b < s.timing_hist.size(); ++b ) {
                if( s.timing_hist[b] == 0 ) continue;
                // Use exact double bounds so the Python midpoint (sqrt(lo*hi)) is accurate.
                double const lo_d = (b == 0) ? 0.0 : std::pow( 10.0, (b - 1) / 10.0 );
                double const hi_d = std::pow( 10.0, b / 10.0 );
                csv_grid_prefix( f, i, grid[i] );
                f << name << "," << lo_d << "," << hi_d << "," << s.timing_hist[b] << "\n";
            }
        }
    }
}

inline void write_offsets_csv(
    std::string const& path,
    std::vector<MutateParams> const& grid,
    std::vector<std::map<std::string, BenchmarkStats>> const& stats
) {
    std::ofstream f( path );
    f << "grid_idx,sub_rate,indel_rate,indel_mean_len,damage_rate,decay_lambda,"
         "aligner,offset_bp,start_count,end_count\n";
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( auto const& [name, s] : stats[i] ) {
            // Merge keys from both offset histograms
            std::set<int32_t> keys;
            for( auto const& [k, _] : s.start_offset_hist ) keys.insert( k );
            for( auto const& [k, _] : s.end_offset_hist   ) keys.insert( k );
            for( int32_t k : keys ) {
                auto const si = s.start_offset_hist.find( k );
                auto const ei = s.end_offset_hist.find( k );
                size_t const sc = ( si != s.start_offset_hist.end() ) ? si->second : 0;
                size_t const ec = ( ei != s.end_offset_hist.end()   ) ? ei->second : 0;
                csv_grid_prefix( f, i, grid[i] );
                f << name << "," << k << "," << sc << "," << ec << "\n";
            }
        }
    }
}

inline void write_scores_csv(
    std::string const& path,
    std::vector<MutateParams> const& grid,
    std::vector<std::map<std::string, BenchmarkStats>> const& stats
) {
    std::ofstream f( path );
    f << "grid_idx,sub_rate,indel_rate,indel_mean_len,damage_rate,decay_lambda,"
         "aligner,score,count\n";
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( auto const& [name, s] : stats[i] ) {
            for( auto const& [score, count] : s.score_hist ) {
                csv_grid_prefix( f, i, grid[i] );
                f << name << "," << score << "," << count << "\n";
            }
        }
    }
}

inline void write_summary_csv(
    std::string const& path,
    std::vector<MutateParams> const& grid,
    std::vector<std::map<std::string, BenchmarkStats>> const& stats
) {
    std::ofstream f( path );
    f << "grid_idx,sub_rate,indel_rate,indel_mean_len,damage_rate,decay_lambda,"
         "aligner,successful,failed,total_ns,mean_ns,min_ns,max_ns\n";
    for( size_t i = 0; i < grid.size(); ++i ) {
        for( auto const& [name, s] : stats[i] ) {
            size_t const   total   = s.successful_alignments + s.failed_alignments;
            uint64_t const mean_ns = total > 0 ? s.total_ns / total : 0;
            uint64_t const min_ns  = ( s.min_ns == std::numeric_limits<uint64_t>::max() ) ? 0 : s.min_ns;
            csv_grid_prefix( f, i, grid[i] );
            f << name << ","
              << s.successful_alignments << "," << s.failed_alignments << ","
              << s.total_ns << "," << mean_ns << "," << min_ns << "," << s.max_ns << "\n";
        }
    }
}

inline void write_windows_csv(
    std::string const& path,
    std::map<size_t, size_t> const& window_hist
) {
    std::ofstream f( path );
    f << "window_size_bp,count\n";
    for( auto const& [size, count] : window_hist ) {
        f << size << "," << count << "\n";
    }
}

inline void write_all_csvs(
    std::string const& bench_dir,
    std::vector<MutateParams> const& grid,
    std::vector<std::vector<size_t>> const& length_hists,
    std::vector<std::map<std::string, BenchmarkStats>> const& stats,
    std::map<size_t, size_t> const& window_hist
) {
    std::filesystem::create_directories( bench_dir );
    write_lengths_csv( bench_dir + "/lengths.csv", grid, length_hists );
    write_timing_csv(  bench_dir + "/timing.csv",  grid, stats        );
    write_offsets_csv( bench_dir + "/offsets.csv", grid, stats        );
    write_scores_csv(  bench_dir + "/scores.csv",  grid, stats        );
    write_summary_csv( bench_dir + "/summary.csv", grid, stats        );
    write_windows_csv( bench_dir + "/windows.csv", window_hist        );
    std::cout << "\nCSV output written to: " << bench_dir << "/\n";
}
