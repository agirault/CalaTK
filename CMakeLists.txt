#---------------------------------------------
# CMake file for the CALATK project STILL CHANGING

cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0012 NEW)
endif(COMMAND cmake_policy)

project(CALATK)

# Version setup
set(CALATK_VERSION_MAJOR "0")
set(CALATK_VERSION_MINOR "1")

mark_as_advanced(Verbose_Compilation)
option(Verbose_Compilation OFF)
if(Verbose_Compilation)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif(Verbose_Compilation)

#
# Output paths
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

#
# What to compile
#
option(COMPILE_APPLICATIONS "Compile command line applications." ON)
option(BUILD_DOCUMENTATION  "Build the documentation." OFF)

# option to use NLOPT as an optimizer
option(USE_NLOPT "Use NLOpt (not recommended)" OFF)
if(USE_NLOPT)
  # So that it finds nlopt
  find_library(NLOPT_LIB nlopt)
  add_definitions(-DUSE_NLOPT)
else(USE_NLOPT)
  set(NLOPT_LIB "")
endif(USE_NLOPT)

# option to use LBFGS as an optimizer
option(USE_LBFGS "Use LBFGS" OFF)
if(USE_LBFGS)
  add_definitions(-DUSE_LBFGS)
endif(USE_LBFGS)

# option to use IPOPT as an optimizer
option(USE_IPOPT "Use IpOpt" OFF)
if(USE_IPOPT)
  # Find the Ipopt library and includes
  find_path(IPOPT_PATH IpoptConfig.h DOC "where IpoptConfig.h is found")
  include_directories(${IPOPT_PATH})
  find_path(IPOPT_LIB_PATH ipopt DOC "path where libipopt is located")

  set(IPOPT_LIBS
    ${IPOPT_LIB_PATH}/libcoinlapack.a
    ${IPOPT_LIB_PATH}/libcoinblas.a
    ${IPOPT_LIB_PATH}/libcoinasl.a
    ${IPOPT_LIB_PATH}/libipopt.a
    ${IPOPT_LIB_PATH}/libipoptamplinterface.a
    ${IPOPT_LIB_PATH}/libcoinmumps.a
    -lgfortran
  )
  add_definitions(-DUSE_IPOPT)
else(USE_IPOPT)
  set(IPOPT_LIBS "")
endif(USE_IPOPT)


#
# ITK Setup
#
option( USE_SYSTEM_ITK "Exclude ITK from SuperBuild and use an existing build instead." OFF )
mark_as_advanced( USE_SYSTEM_ITK )
if( USE_SYSTEM_ITK )
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
  find_path(FFTW_PATH fftw3.h DOC "include directory for FFTW")
  include_directories(${FFTW_PATH})
  find_library(FFTW_LIB fftw3)
  find_library(FFTWF_LIB fftw3f)
endif()

# Superbuild setup
option( CALATK_USE_SUPERBUILD
  "Build calatk and the projects it depends on via SuperBuild.cmake." ON )
mark_as_advanced( CALATK_USE_SUPERBUILD )
if( CALATK_USE_SUPERBUILD AND WIN32 )
  # FFTW is not built.
  message( WARNING "Superbuild is not supported on Windows yet." )
endif()

if( CALATK_USE_SUPERBUILD AND NOT WIN32 )
  include( "${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake" )
  return()
else()
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
  find_path(FFTW_PATH fftw3.h DOC "include directory for FFTW")
  include_directories(${FFTW_PATH})
  find_library(FFTW_LIB fftw3)
  find_library(FFTWF_LIB fftw3f)
endif()

#
# Shared libraries do not link correctly yet on Windows.
#
if(WIN32)
  set(BUILD_SHARED_LIBS OFF)
endif()
  set(BUILD_SHARED_LIBS OFF)
#
# Configure calatkConfigure.h.in
#
set(CALATK_SYSTEM_UNIX ${UNIX})
set(CALATK_SYSTEM_WINDOWS ${WIN32})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions(-DDARWIN)
endif()

configure_file(${CMAKE_SOURCE_DIR}/calatkConfigure.h.in ${CMAKE_BINARY_DIR}/calatkConfigure.h)

#
# Add OpenMP flag
#
if(WIN32)
  if(CMAKE_GENERATOR MATCHES Visual* )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    set(COMPILER_HAS_OPENMP TRUE)
  #else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
  endif()
else(CMAKE_COMPILER_IS_GNUCC)
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_FLAGS "-fopenmp -Werror")
  check_cxx_source_compiles("int main() { return 0; }" COMPILER_HAS_OPENMP)
  if(COMPILER_HAS_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif()
endif()

# Suppress duplicate dynlib warnings.
# For more information, see http://public.kitware.com/Bug/view.php?id=10179
#
#if(APPLE)
#  set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -w)
#endif()

#
# Allow for boost support for some of the libraries (this is optional)
#
option(USE_BOOST "Use boost." OFF)
if(USE_BOOST)
  find_path(BOOST_INCLUDE_DIR boost/config.hpp)
  find_path(BOOST_LIBRARY_DIR boost_date_time-vc71-mt.lib)

  link_directories(${BOOST_LIBRARY_DIR})
  include_directories(${BOOST_INCLUDE_DIR})

  # Win32 has automatic linking of boost libraries so only add boost
  # libaries on unix
  if(NOT WIN32 OR CYGWIN)
    set(BOOST_LIBRARIES boost_program_options)
  endif(NOT WIN32 OR CYGWIN)
  add_definitions(-DUSE_BOOST)
endif(USE_BOOST)

#
# Allow to choose the floating point type
#
option(CHOOSE_FLOATING_POINT_TYPE "Choose floating point type on command line." ON)
if(CHOOSE_FLOATING_POINT_TYPE)
  add_definitions(-DFLOATING_POINT_CHOICE)
endif(CHOOSE_FLOATING_POINT_TYPE)

#
# Include directories
#
set(CALATK_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Algorithms
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Base
  ${CMAKE_SOURCE_DIR}/Code/Libraries/DataStructures
  ${CMAKE_SOURCE_DIR}/Code/Libraries/ImageManager
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Utilities
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Numerics
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Metrics
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Kernels
  ${CMAKE_SOURCE_DIR}/Code/Libraries/States
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Resamplers
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Interfaces
  ${CMAKE_SOURCE_DIR}/Code/Libraries/Interpolators
  ${CMAKE_SOURCE_DIR}/Code/Libraries/ObjectiveFunctions
  ${CMAKE_SOURCE_DIR}/Code/External/ITKTestingCommands
  ${CMAKE_SOURCE_DIR}/Code/External/jsoncpp
  ${CMAKE_BINARY_DIR}
  # Casey's Code
  ${CMAKE_SOURCE_DIR}/Utilities/DTIProcess/Library
)


if(USE_NLOPT)
  include_directories(${CMAKE_SOURCE_DIR}/Code/External/nlopt-2.2.4/api)
endif(USE_NLOPT)

if(USE_LBFGS)
  include_directories(${CMAKE_SOURCE_DIR}/Code/External/libbfgs)
endif(USE_LBFGS)

include_directories(${CALATK_INCLUDE_DIRS})

#
# Configure CALATKConfig.cmake.in and UseCALATK.cmake.in for external projects
#
set(CALATK_INCLUDE_DIRS_CONFIG ${CALATK_INCLUDE_DIRS} ${ITK_INCLUDE_DIRS})
set(CALATK_LIBRARY_DIRS_CONFIG ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${ITK_LIBRARY_DIRS})
set(CALATK_USE_FILE ${CMAKE_BINARY_DIR}/UseCALATK.cmake)

configure_file(${CMAKE_SOURCE_DIR}/UseCALATK.cmake.in
               ${CALATK_USE_FILE} COPYONLY IMMEDIATE)
configure_file(${CMAKE_SOURCE_DIR}/CALATKConfig.cmake.in
               ${CMAKE_BINARY_DIR}/CALATKConfig.cmake)


#
# Create the documentation
#
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/Utilities/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile)
    add_custom_target(Documentation ALL doxygen ${CMAKE_BINARY_DIR}/Doxyfile)
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCUMENTATION)

# Set up for testing
include(CTest)
enable_testing()

# For Cppcheck.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CALATK_SOURCE_DIR}/Scripts/Cppcheck" )

# Build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSEInstaller.txt")
set(CPACK_RESOURCE_FILE_README
     "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_WELCOME
     "${CMAKE_CURRENT_SOURCE_DIR}/WELCOME.txt")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Image Registration and Atlas-Building")
set(CPACK_PACKAGE_VENDOR "The CALATK Team")

set(CPACK_PACKAGE_VERSION_MAJOR "${CALATK_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${CALATK_VERSION_MINOR}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CALATK - Cross-sectional and longitudinal atlas-building toolkit")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "CALATK-${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

include (CPack)
add_definitions(-DCALATK_VERSION_MAJOR=${CALATK_VERSION_MAJOR} -DCALATK_VERSION_MINOR=${CALATK_VERSION_MINOR})

#
# Go to Code
#
add_subdirectory(${CMAKE_SOURCE_DIR}/Code)

#
# Configure CTestCustom.cmake.in
#
configure_file(${CALATK_SOURCE_DIR}/Scripts/CTestCustom.cmake.in
  ${CALATK_BINARY_DIR}/CTestCustom.cmake)


