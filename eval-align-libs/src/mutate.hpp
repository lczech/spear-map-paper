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
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// =================================================================================================
//     MutateParams
// =================================================================================================

struct MutateParams
{
    double sub_rate       = 0.0;   // per-base substitution probability (models divergence)
    double indel_rate     = 0.0;   // per-base probability of an indel event
    double indel_mean_len = 1.5;   // mean indel length (geometric distribution, min 1)
    double damage_rate    = 0.0;   // C→T / G→A probability at position 0 from each terminus
    double decay_lambda   = 0.3;   // exponential decay rate for damage away from termini
    size_t damage_tail    = 10;    // number of positions from each end to apply damage to
};

// =================================================================================================
//     Helpers
// =================================================================================================

namespace {

std::string reverse_complement_mutated( std::string const& seq )
{
    std::string rc( seq.size(), 'N' );
    for( size_t i = 0; i < seq.size(); ++i ) {
        size_t const rev = seq.size() - 1 - i;
        switch( seq[i] ) {
            case 'A': rc[rev] = 'T'; break;
            case 'T': rc[rev] = 'A'; break;
            case 'C': rc[rev] = 'G'; break;
            case 'G': rc[rev] = 'C'; break;
            default:  rc[rev] = 'N'; break;
        }
    }
    return rc;
}

} // namespace

// =================================================================================================
//     mutate
// =================================================================================================

// Returns a copy of `read` with `forward` mutated and `reverse` regenerated.
// All other fields (origin coords, window) are copied unchanged.
// Order: divergence (substitutions then indels) first, then aDNA damage.
ShatteredRead mutate(
    ShatteredRead const& read,
    MutateParams const& params,
    std::mt19937_64& rng
) {
    ShatteredRead result = read;
    std::string& seq = result.forward;
    size_t len = seq.size();

    constexpr char bases[4] = { 'A', 'C', 'G', 'T' };

    std::uniform_int_distribution<int>  base4( 0, 3 );  // random base for insertions
    std::uniform_int_distribution<int>  base3( 0, 2 );  // offset into 3 non-current bases
    std::uniform_real_distribution<>    prob( 0.0, 1.0 );

    // --- Step 1: Substitutions ---
    // Draw count from Binomial(len, sub_rate), then sample positions without replacement.

    if( params.sub_rate > 0.0 && len > 0 ) {
        std::binomial_distribution<int> binom(
            static_cast<int>( len ), params.sub_rate
        );
        int const k = binom( rng );
        if( k > 0 ) {
            // Sample k unique positions using std::sample
            std::vector<size_t> indices( len );
            std::iota( indices.begin(), indices.end(), size_t{0} );
            std::vector<size_t> positions( static_cast<size_t>( k ) );
            std::sample( indices.begin(), indices.end(), positions.begin(), k, rng );

            for( size_t pos : positions ) {
                char& c = seq[pos];
                int cur = 0;
                for( ; cur < 4; ++cur ) {
                    if( bases[cur] == c ) break;
                }
                c = bases[ (cur + 1 + base3( rng )) % 4 ];
            }
        }
    }

    // --- Step 2: Indels ---
    // Draw count from Binomial(len, indel_rate), sample positions, then apply in
    // reverse order (deletions) or accumulate into a new string (both types).

    if( params.indel_rate > 0.0 && len > 0 ) {
        std::binomial_distribution<int> binom(
            static_cast<int>( len ), params.indel_rate
        );
        int const k = binom( rng );
        if( k > 0 ) {
            std::geometric_distribution<int> geo(
                params.indel_mean_len > 0.0 ? 1.0 / params.indel_mean_len : 1.0
            );

            // Sample k positions and sort them so we can process left to right.
            std::vector<size_t> indices( len );
            std::iota( indices.begin(), indices.end(), size_t{0} );
            std::vector<size_t> positions( static_cast<size_t>( k ) );
            std::sample( indices.begin(), indices.end(), positions.begin(), k, rng );
            std::sort( positions.begin(), positions.end() );

            std::string mutated;
            mutated.reserve( seq.size() );

            size_t prev = 0;
            for( size_t pos : positions ) {
                // Copy unchanged bases up to this event position
                mutated.append( seq, prev, pos - prev );

                int const ilen = 1 + geo( rng );
                if( prob( rng ) < 0.5 ) {
                    // Insertion: random bases inserted before seq[pos]
                    for( int j = 0; j < ilen; ++j ) {
                        mutated += bases[ base4( rng ) ];
                    }
                    mutated += seq[pos];   // keep the base at this position
                } else {
                    // Deletion: skip ilen bases starting at pos
                    pos += static_cast<size_t>( ilen ) - 1;  // -1 because prev advances past pos
                }
                prev = pos + 1;
            }
            // Copy remainder
            if( prev < seq.size() ) {
                mutated.append( seq, prev, seq.size() - prev );
            }

            seq = std::move( mutated );
            len = seq.size();
        }
    }

    // --- Step 3: aDNA damage ---
    // Only process the first and last `damage_tail` positions — the exponential
    // decay makes probabilities negligible beyond that range.

    if( params.damage_rate > 0.0 && len > 0 ) {
        size_t const tail = std::min( params.damage_tail, len / 2 );

        // 5' end: C→T
        for( size_t j = 0; j < tail; ++j ) {
            if( seq[j] == 'C' ) {
                double const p = params.damage_rate
                    * std::exp( -static_cast<double>( j ) * params.decay_lambda );
                if( prob( rng ) < p ) {
                    seq[j] = 'T';
                }
            }
        }
        // 3' end: G→A
        for( size_t j = 0; j < tail; ++j ) {
            size_t const pos = len - 1 - j;
            if( seq[pos] == 'G' ) {
                double const p = params.damage_rate
                    * std::exp( -static_cast<double>( j ) * params.decay_lambda );
                if( prob( rng ) < p ) {
                    seq[pos] = 'A';
                }
            }
        }
    }

    result.reverse = reverse_complement_mutated( seq );
    return result;
}
