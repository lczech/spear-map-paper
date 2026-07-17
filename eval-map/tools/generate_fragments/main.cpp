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

#include "mutate.hpp"
#include "shatter.hpp"

#include "genesis/genesis.hpp"
#include <CLI/CLI.hpp>

#include <cstdint>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>

using namespace genesis::sequence;
using namespace genesis::util;
using namespace genesis::util::io;

// ================================================================================================
//     main
// ================================================================================================

int main( int argc, char** argv )
{
    // Activate logging
    genesis::util::core::Logging::log_to_stdout();
    genesis::util::core::Logging::details.date = true;
    genesis::util::core::Logging::details.time = true;

    // Genesis' gzip output stream requires the global thread pool to be initialized, even though
    // this tool's own logic (shatter/mutate/write) is purely sequential.
    genesis::util::core::Options::get().init_global_thread_pool( 1 );

    // -------------------------------------------------------------------
    //     cli setup
    // -------------------------------------------------------------------

    CLI::App app{
        "Simulate reads (shattered and mutated fragments) from a reference genome, "
        "paired with their known ground-truth origin, for use in mapping benchmarks."
    };

    std::string reference;
    std::string out;
    std::string truth_out;
    unsigned int seed = 42;

    ShatterParams shatter_params;
    MutateParams  mutate_params;

    app.add_option( "--reference", reference, "Reference genome FASTA file" )
       ->required();
    app.add_option( "--out", out, "Output FASTQ file (gzipped)" )
       ->required();
    app.add_option( "--truth-out", truth_out,
        "Output truth TSV file (columns: read_id, chrom, start, end, strand)" )
       ->required();

    app.add_option( "--min-len", shatter_params.min_len, "Minimum fragment length" )
       ->capture_default_str();
    app.add_option( "--max-len", shatter_params.max_len, "Maximum fragment length" )
       ->capture_default_str();
    app.add_option( "--gamma-shape", shatter_params.gamma_a, "Gamma distribution shape" )
       ->capture_default_str();
    app.add_option( "--gamma-scale", shatter_params.gamma_b, "Gamma distribution scale" )
       ->capture_default_str();

    app.add_option( "--sub-rate", mutate_params.sub_rate,
        "Per-base substitution probability, modeling divergence from the reference" )
       ->capture_default_str()
       ->check( CLI::Range( 0.0, 1.0 ));
    app.add_option( "--indel-rate", mutate_params.indel_rate,
        "Per-base probability of an indel event" )
       ->capture_default_str()
       ->check( CLI::Range( 0.0, 1.0 ));
    app.add_option( "--indel-mean-len", mutate_params.indel_mean_len,
        "Mean indel length (geometric distribution, min 1)" )
       ->capture_default_str();
    app.add_option( "--damage-rate", mutate_params.damage_rate,
        "aDNA damage probability at read termini (C->T at 5', G->A at 3')" )
       ->capture_default_str()
       ->check( CLI::Range( 0.0, 1.0 ));
    app.add_option( "--decay-lambda", mutate_params.decay_lambda,
        "Exponential decay rate of damage probability away from the termini" )
       ->capture_default_str();
    app.add_option( "--damage-tail", mutate_params.damage_tail,
        "Number of positions from each end subject to damage" )
       ->capture_default_str();

    app.add_option( "--seed", seed, "Random seed" )
       ->capture_default_str();

    // Parse
    CLI11_PARSE( app, argc, argv );
    shatter_params.seed = seed;

    LOG_MSG << "Reference:      " << reference << "\n"
            << "Min len:        " << shatter_params.min_len << "\n"
            << "Max len:        " << shatter_params.max_len << "\n"
            << "Gamma shape:    " << shatter_params.gamma_a << "\n"
            << "Gamma scale:    " << shatter_params.gamma_b << "\n"
            << "Sub rate:       " << mutate_params.sub_rate << "\n"
            << "Indel rate:     " << mutate_params.indel_rate << "\n"
            << "Indel mean len: " << mutate_params.indel_mean_len << "\n"
            << "Damage rate:    " << mutate_params.damage_rate << "\n"
            << "Decay lambda:   " << mutate_params.decay_lambda << "\n"
            << "Damage tail:    " << mutate_params.damage_tail << "\n"
            << "Seed:           " << seed << "\n"
            << "Output:         " << out << "\n"
            << "Truth output:   " << truth_out << std::endl;

    // Randomness. Reused for both the mutation model and the strand coin flip below,
    // so that the whole run is reproducible from a single --seed.
    std::mt19937_64 rng( seed );
    std::uniform_real_distribution<double> strand_coin( 0.0, 1.0 );

    // -------------------------------------------------------------------
    //     output prep
    // -------------------------------------------------------------------

    FastqWriter writer;
    writer.fill_missing_quality( 40 );
    auto fastq_stream = FastxOutputStream( io::to_gzip_file( out ), writer );

    std::ofstream truth_stream( truth_out );
    if( !truth_stream ) {
        throw std::runtime_error( "Could not open truth output file: " + truth_out );
    }
    truth_stream << "read_id\tchrom\tstart\tend\tstrand\n";

    // -------------------------------------------------------------------
    //     main loop
    // -------------------------------------------------------------------

    size_t num_reads    = 0;
    size_t num_filtered = 0;

    for( auto const& read : shatter( reference, shatter_params )) {
        // Fragments touching assembly gaps have no well-defined true origin to score against.
        if( !passes_filter( read )) {
            ++num_filtered;
            continue;
        }

        auto const mutated = mutate( read, mutate_params, rng );

        // Randomly emit the fragment from either genomic strand, as real sequencing reads do,
        // so that aligners' reverse-complement mapping is exercised as well.
        bool const is_reverse = strand_coin( rng ) < 0.5;
        std::string const& seq = is_reverse ? mutated.reverse : mutated.forward;
        char const strand = is_reverse ? '-' : '+';

        auto const read_id =
            mutated.chromosome + "_" + std::to_string( mutated.origin_start )
        ;

        fastq_stream << Sequence( read_id, seq );
        truth_stream
            << read_id << '\t'
            << mutated.chromosome << '\t'
            << mutated.origin_start << '\t'
            << mutated.origin_end << '\t'
            << strand << '\n'
        ;

        ++num_reads;
    }

    // DONE!
    LOG_MSG << "num_reads: " << num_reads;
    LOG_MSG << "num_filtered (N-containing): " << num_filtered;
    LOG_MSG << "finished";
    return 0;
}
