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
    set( ColorYellow "${Esc}[33m" )
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

project( eval-align-libs CXX C )

# --------------------------------------------------------------------------------------------------
#   Build Type
# --------------------------------------------------------------------------------------------------

if( CMAKE_BUILD_TYPE STREQUAL "" )
    set( CMAKE_BUILD_TYPE RELEASE )
endif()
message( STATUS "${ColorBlue}eval-align-libs build type: ${CMAKE_BUILD_TYPE}${ColorEnd}" )

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

set( DEPS_DIR "${PROJECT_SOURCE_DIR}/../deps" )

# --- Genesis ------------------------------------------------------------------

set( GENESIS_CPP_STANDARD        "20" CACHE BOOL "" FORCE )
set( GENESIS_BUILD_APPLICATIONS  OFF  CACHE BOOL "" FORCE )
set( GENESIS_BUILD_PYTHON_MODULE OFF  CACHE BOOL "" FORCE )
set( GENESIS_BUILD_TESTS         OFF  CACHE BOOL "" FORCE )
set( GENESIS_USE_HTSLIB          OFF  CACHE BOOL "" FORCE )

set( genesis_SUBMODULE_DIR "${DEPS_DIR}/genesis" )
message( STATUS "Looking for genesis" )

if( EXISTS "${genesis_SUBMODULE_DIR}/CMakeLists.txt" )
    message( STATUS "${ColorGreen}Found genesis submodule at ${genesis_SUBMODULE_DIR}${ColorEnd}" )
    add_subdirectory( "${genesis_SUBMODULE_DIR}" "${CMAKE_BINARY_DIR}/genesis" )
    set( genesis_SOURCE "${genesis_SUBMODULE_DIR}" )
else()
    message( FATAL_ERROR
        "${ColorRed}Genesis submodule not found at ${genesis_SUBMODULE_DIR}. "
        "Run: git submodule update --init --recursive${ColorEnd}"
    )
endif()

# Genesis headers — marked SYSTEM to suppress warnings from their headers
include_directories( SYSTEM "${genesis_SOURCE}/lib" )

# --------------------------------------------------------------------------------------------------
#   Target
# --------------------------------------------------------------------------------------------------

add_executable( eval-align-libs src/main.cpp )

target_include_directories( eval-align-libs PRIVATE "${PROJECT_SOURCE_DIR}/src" )

target_link_libraries( eval-align-libs PRIVATE genesis_lib_static )

if( GENESIS_HAS_LTO )
    message( STATUS "Building with LTO/IPO support" )
    set_property( TARGET eval-align-libs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE )
else()
    message( STATUS "Building without LTO/IPO support" )
endif()
