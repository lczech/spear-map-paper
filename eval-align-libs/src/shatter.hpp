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

#include "generator.hpp"

#include "genesis/genesis.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>

// =================================================================================================
//     ShatteredRead
// =================================================================================================

// A fragment extracted from a reference genome, along with the reference window it will be
// aligned against in the benchmark.
struct ShatteredRead
{
    std::string chromosome;
    size_t      origin_start = 0;   // 0-based, inclusive
    size_t      origin_end   = 0;   // 0-based, exclusive

    std::string forward;            // read sequence (uppercase; non-ACGT bases set to N)
    std::string reverse;            // reverse complement of forward

    size_t      window_start = 0;   // origin_start rounded down to interval boundary
    size_t      window_end   = 0;   // origin_end   rounded up   to interval boundary
    std::string window;             // reference sequence in [window_start, window_end)
};

// =================================================================================================
//     ShatterParams
// =================================================================================================

struct ShatterParams
{
    double   gamma_a        = 2.5;   // gamma shape  \  together give a peak around 55-65 bp,
    double   gamma_b        = 30.0;  // gamma scale  /  range [min_len, max_len]
    int      min_len        = 30;
    int      max_len        = 150;
    size_t   interval_size  = 128;   // reference window is aligned to multiples of this
    size_t   window_padding = 16;    // extra reference context added outside the interval window
    uint64_t seed           = 42;
};

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
//     Helpers
// =================================================================================================

namespace {

// Uppercases seq in place, replaces non-ACGT with N, and returns the reverse complement.
std::string clean_and_rc( std::string& seq )
{
    std::string rc( seq.size(), 'N' );
    for( size_t i = 0; i < seq.size(); ++i ) {
        seq[i] = static_cast<char>( std::toupper( static_cast<unsigned char>( seq[i] )));
        size_t const rev = seq.size() - 1 - i;
        switch( seq[i] ) {
            case 'A': rc[rev] = 'T'; break;
            case 'T': rc[rev] = 'A'; break;
            case 'C': rc[rev] = 'G'; break;
            case 'G': rc[rev] = 'C'; break;
            default:  seq[i] = 'N'; rc[rev] = 'N'; break;
        }
    }
    return rc;
}

} // namespace

// =================================================================================================
//     shatter
// =================================================================================================

// Lazily yields ShatteredRead fragments from a FASTA file, one at a time.
// Read lengths follow a gamma distribution clamped to [min_len, max_len].
// Each read is paired with an interval-aligned reference window.
// Fragments whose window would extend past the chromosome end are skipped.
// N-filtering is left to the caller (see harness.hpp).
Generator<ShatteredRead> shatter(
    std::string const& fasta_path,
    ShatterParams const& params = ShatterParams{}
) {
    using namespace genesis::sequence;
    using namespace genesis::util;
    using namespace genesis::util::io;

    std::mt19937_64 rng( params.seed );
    std::gamma_distribution<double> gamma_dist( params.gamma_a, params.gamma_b );

    LOG_MSG << "Shattering " << fasta_path;

    auto fasta_stream = FastaInputStream( from_file( fasta_path ));
    for( auto const& seq : fasta_stream ) {
        auto const chr      = seq.label().substr( 0, seq.label().find( ' ' ));
        size_t const chr_len = seq.size();
        LOG_MSG1 << "Processing chromosome " << chr << " (length " << chr_len << ")";

        size_t pos = 0;
        while( pos < chr_len ) {

            // Draw a length in [min_len, max_len] from the gamma distribution.
            int len = 0;
            do {
                len = static_cast<int>( std::round( gamma_dist( rng )));
            } while( len < params.min_len || len > params.max_len );

            // Clamp to remaining chromosome; break if too short after clamping.
            if( pos + static_cast<size_t>( len ) > chr_len ) {
                len = static_cast<int>( chr_len - pos );
            }
            if( len < params.min_len ) {
                break;
            }

            // Compute interval-aligned reference window boundaries.
            size_t const win_start_ivl =
                ( pos / params.interval_size ) * params.interval_size;
            size_t const win_end_ivl =
                (( pos + len + params.interval_size - 1 ) / params.interval_size )
                * params.interval_size;

            // Skip fragments whose interval window would extend past the chromosome end.
            if( win_end_ivl > chr_len ) {
                pos += len;
                continue;
            }

            // Apply padding, clamped to chromosome boundaries.
            size_t const win_start = win_start_ivl > params.window_padding
                ? win_start_ivl - params.window_padding : 0;
            size_t const win_end = std::min( chr_len, win_end_ivl + params.window_padding );

            // Extract and clean the read sequence.
            std::string fwd = seq.sites().substr( pos, len );
            std::string rev = clean_and_rc( fwd );

            // Extract and uppercase the reference window.
            std::string win = seq.sites().substr( win_start, win_end - win_start );
            genesis::util::text::to_upper_inplace( win );

            co_yield ShatteredRead{
                chr,
                pos,
                pos + static_cast<size_t>( len ),
                std::move( fwd ),
                std::move( rev ),
                win_start,
                win_end,
                std::move( win )
            };

            pos += len;
        }
    }

    LOG_MSG << "Finished shattering " << fasta_path;
}
