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

#include "stats.hpp"

#include "edlib.h"

#include <string>

// =================================================================================================
//     edlib aligner
// =================================================================================================

// Semi-global (infix / HW) alignment of query against target using edlib.
// No precomputation step exists for edlib, so there is no hot/cold split —
// this single function covers the full per-call cost.
// Returns edit distance as score; start and end are 0-based positions in target.
inline AlignResult align_edlib( std::string const& query, std::string const& target )
{
    EdlibAlignResult raw = edlibAlign(
        query.c_str(),  static_cast<int>( query.size() ),
        target.c_str(), static_cast<int>( target.size() ),
        edlibNewAlignConfig( -1, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0 )
    );

    AlignResult result;
    if( raw.status != EDLIB_STATUS_OK || raw.editDistance == -1 || raw.numLocations == 0 ) {
        result.failed = true;
    } else {
        result.start = raw.startLocations[0];
        result.end   = raw.endLocations[0] + 1;   // edlib end is inclusive; convert to exclusive
        result.score = raw.editDistance;
    }

    edlibFreeAlignResult( raw );
    return result;
}
