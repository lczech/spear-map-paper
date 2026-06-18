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

#include "harness.hpp"
#include "shatter.hpp"

#include "genesis/util/core/logging.hpp"

#include <iostream>

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
    //     Generate reads
    // -------------------------------------------------------------------------

    LOG_MSG << "Alignment library benchmarking";

    auto const stats = collect_stats( shatter( fasta_path ));
    print_stats( stats );

    return 0;
}
