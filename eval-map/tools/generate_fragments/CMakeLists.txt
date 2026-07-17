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

set( CMAKE_DISABLE_IN_SOURCE_BUILD ON )
set( CMAKE_DISABLE_SOURCE_CHANGES  ON )

if( NOT ( CMAKE_VERSION VERSION_LESS 3.0 ))
    cmake_policy( SET CMP0025 NEW )
    cmake_policy( SET CMP0028 NEW )
endif()

# --------------------------------------------------------------------------------------------------
#   Project Definition
# --------------------------------------------------------------------------------------------------

project( generate_fragments CXX C )

if( CMAKE_BUILD_TYPE STREQUAL "" )
    set( CMAKE_BUILD_TYPE RELEASE )
endif()
message( STATUS "${ColorBlue}generate_fragments build type: ${CMAKE_BUILD_TYPE}${ColorEnd}" )

set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

if( MSVC OR ( CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" ))
    set( WARN_FLAGS "-Wall" )
else()
    set( WARN_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors" )
    if( CMAKE_BUILD_TYPE MATCHES DEBUG )
        set( WARN_FLAGS "${WARN_FLAGS} -Werror" )
    endif()
endif()

set( CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++20 ${WARN_FLAGS}" )
set( CMAKE_CXX_FLAGS_DEBUG   "-O2 -DDEBUG -g -ggdb3 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
set( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" )

if( APPLE )
    set( CMAKE_MACOSX_RPATH ON )
endif()

if( NOT ( CMAKE_VERSION VERSION_LESS 3.9 ))
    cmake_policy( SET CMP0069 NEW )
endif()

# --------------------------------------------------------------------------------------------------
#   Dependencies
# --------------------------------------------------------------------------------------------------

# --- eval_common (genesis, shatter/mutate read simulation) ---------------------

set( eval_common_DIR "${PROJECT_SOURCE_DIR}/../../../common" )
message( STATUS "Looking for common eval tooling" )

if( EXISTS "${eval_common_DIR}/CMakeLists.txt" )
    message( STATUS "${ColorGreen}Found common eval tooling at ${eval_common_DIR}${ColorEnd}" )
    add_subdirectory( "${eval_common_DIR}" "${CMAKE_BINARY_DIR}/common" )
else()
    message( FATAL_ERROR
        "${ColorRed}common eval tooling not found at ${eval_common_DIR}.${ColorEnd}"
    )
endif()

# --------------------------------------------------------------------------------------------------
#   Target
# --------------------------------------------------------------------------------------------------

add_executable( generate_fragments main.cpp )

target_link_libraries( generate_fragments PRIVATE eval_common )

if( GENESIS_HAS_LTO )
    message( STATUS "Building with LTO/IPO support" )
    set_property( TARGET generate_fragments PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE )
else()
    message( STATUS "Building without LTO/IPO support" )
endif()
