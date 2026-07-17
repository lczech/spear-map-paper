# SPEAR - Sorting Petabytes of Environmental and Ancient Reads
# Copyright (C) 2026 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@sund.ku.dk>
# University of Copenhagen, Globe Institute, Section for GeoGenetics
# Oster Voldgade 5-7, 1350 Copenhagen K, Denmark

# --------------------------------------------------------------------------------------------------
#   CMake Init
# --------------------------------------------------------------------------------------------------

cmake_minimum_required( VERSION 3.14 FATAL_ERROR )

if( NOT WIN32 )
    string( ASCII 27 Esc )
    set( ColorBlue   "${Esc}[34m" )
    set( ColorGreen  "${Esc}[32m" )
    set( ColorRed    "${Esc}[31m" )
    set( ColorEnd    "${Esc}[0m"  )
endif()

# --------------------------------------------------------------------------------------------------
#   Genesis
# --------------------------------------------------------------------------------------------------

# This is the single place where the genesis submodule is located and configured.
# Consumers of the `eval_common` target below get genesis (and its bundled CLI11) transitively,
# and do not need to look for or add_subdirectory() genesis themselves.

set( GENESIS_CPP_STANDARD        "20" CACHE BOOL "" FORCE )
set( GENESIS_BUILD_APPLICATIONS  OFF  CACHE BOOL "" FORCE )
set( GENESIS_BUILD_PYTHON_MODULE OFF  CACHE BOOL "" FORCE )
set( GENESIS_BUILD_TESTS         OFF  CACHE BOOL "" FORCE )
set( GENESIS_USE_HTSLIB          OFF  CACHE BOOL "" FORCE )

set( genesis_SUBMODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../deps/genesis" )
message( STATUS "Looking for genesis" )

if( EXISTS "${genesis_SUBMODULE_DIR}/CMakeLists.txt" )
    message( STATUS "${ColorGreen}Found genesis submodule at ${genesis_SUBMODULE_DIR}${ColorEnd}" )
    add_subdirectory( "${genesis_SUBMODULE_DIR}" "${CMAKE_BINARY_DIR}/genesis" )

    # Genesis exports GENESIS_HAS_LTO to its parent scope (this file). Bubble it one more level
    # up, so that consumers of eval_common (which add_subdirectory() this file) can still check
    # it to decide whether to enable IPO/LTO on their own targets.
    set( GENESIS_HAS_LTO ${GENESIS_HAS_LTO} PARENT_SCOPE )
else()
    message( FATAL_ERROR
        "${ColorRed}Genesis submodule not found at ${genesis_SUBMODULE_DIR}. "
        "Run: git submodule update --init --recursive${ColorEnd}"
    )
endif()

# --------------------------------------------------------------------------------------------------
#   eval_common target
# --------------------------------------------------------------------------------------------------

# Header-only shared toolkit for the eval-* benchmarks in this repo (read simulation via
# shatter()/mutate(), the Generator<T> coroutine helper, etc). Depends only on genesis.

add_library( eval_common INTERFACE )

target_include_directories( eval_common INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Genesis headers — marked SYSTEM to suppress warnings from their headers.
target_include_directories( eval_common SYSTEM INTERFACE
    "${genesis_SUBMODULE_DIR}/lib"
)

target_link_libraries( eval_common INTERFACE genesis_lib_static )
