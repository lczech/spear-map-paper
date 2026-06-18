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

#include <cmath>
#include <cstddef>
#include <random>
#include <string>

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

    constexpr char bases[4] = { 'A', 'C', 'G', 'T' };

    std::uniform_real_distribution<double> prob( 0.0, 1.0 );
    std::uniform_int_distribution<int>     base4( 0, 3 );   // random base for insertions
    std::uniform_int_distribution<int>     base3( 0, 2 );   // offset into 3 non-current bases
    std::geometric_distribution<int>       geo(
        params.indel_mean_len > 0.0 ? 1.0 / params.indel_mean_len : 1.0
    );

    // --- Step 1: Substitutions ---

    if( params.sub_rate > 0.0 ) {
        for( char& c : seq ) {
            if( prob( rng ) < params.sub_rate ) {
                int cur = 0;
                for( ; cur < 4; ++cur ) {
                    if( bases[cur] == c ) break;
                }
                c = bases[ (cur + 1 + base3( rng )) % 4 ];
            }
        }
    }

    // --- Step 2: Indels ---
    // Process each original position; insertions add random bases before it,
    // deletions skip ahead. Operates on the already-substituted sequence.

    if( params.indel_rate > 0.0 ) {
        std::string mutated;
        mutated.reserve( seq.size() );

        for( size_t i = 0; i < seq.size(); ) {
            if( prob( rng ) < params.indel_rate ) {
                int const len = 1 + geo( rng );
                if( prob( rng ) < 0.5 ) {
                    // Insertion: insert `len` random bases before current position
                    for( int k = 0; k < len; ++k ) {
                        mutated += bases[ base4( rng ) ];
                    }
                    // fall through to also emit seq[i]
                } else {
                    // Deletion: skip `len` bases
                    i += static_cast<size_t>( len );
                    continue;
                }
            }
            if( i < seq.size() ) {
                mutated += seq[i++];
            }
        }
        seq = std::move( mutated );
    }

    // --- Step 3: aDNA damage ---
    // C→T decaying exponentially from 5' end; G→A from 3' end.

    if( params.damage_rate > 0.0 ) {
        size_t const len = seq.size();
        for( size_t j = 0; j < len; ++j ) {
            double const p5 = params.damage_rate
                * std::exp( -static_cast<double>( j ) * params.decay_lambda );
            double const p3 = params.damage_rate
                * std::exp( -static_cast<double>( len - 1 - j ) * params.decay_lambda );
            if( seq[j] == 'C' && prob( rng ) < p5 ) {
                seq[j] = 'T';
            } else if( seq[j] == 'G' && prob( rng ) < p3 ) {
                seq[j] = 'A';
            }
        }
    }

    result.reverse = reverse_complement_mutated( seq );
    return result;
}
