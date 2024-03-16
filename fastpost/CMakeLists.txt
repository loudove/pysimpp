cmake_minimum_required (VERSION 3.5)

project (fastpost)
enable_language (Fortran)

# make sure that the default is RELEASE build
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE 
      STRING "Choose the type of build, options are: Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Set traceback in fortran compiler 
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -traceback")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbacktrace")
endif()

# The version number
set (FASTPOST_VERSION_MAJOR 1)
set (FASTPOST_VERSION_MINOR 5)
 
# include direcories
include_directories ( ${CMAKE_BINARY_DIR} )
include_directories ( ${PROJECT_BINARY_DIR} )

# lapack sources
file ( GLOB LAPACKSOURCES "src/lapack/*.f")

# core sources 
set ( CORESOURCES
    src/core/argsorting_mod.F90
    src/core/vector3d_mod.F90
    src/core/domain3d_mod.F90
    src/core/lj_mod.F90 )

# interface source
set ( LIBSOURCES 
    src/fastpost.f90 
    src/fastbonded.f90)

# message ( STATUS "PROJECT_SOURCE_DIR (fastpost) : ${PROJECT_SOURCE_DIR}" )
# message ( STATUS "PROJECT_BINARY_DIR (fastpost) : ${PROJECT_BINARY_DIR}" )
# message ( STATUS "SYSTEM NAME (fastpost) : ${CMAKE_SYSTEM_NAME}")

# add core library
add_library ( fastcore STATIC ${LAPACKSOURCES} ${CORESOURCES} )
set_property ( TARGET fastcore PROPERTY POSITION_INDEPENDENT_CODE ON)

# build python module (use -DPYTHON_EXECUTABLE to force a specific interpreter)
# find_package ( PythonInterp 3.6 REQUIRED)
# find_package ( PythonLibs 3.6 REQUIRED)
find_package(Python3 COMPONENTS Interpreter)

# check also: 
# https://github.com/scikit-build/scikit-build/blob/master/skbuild/resources/cmake/FindF2PY.cmake
#
find_program ( F2PY_EXECUTABLE NAMES f2py${PYTHON_VERSION_MAJOR} f2py REQUIRED)
message ( STATUS "Found f2py: ${F2PY_EXECUTABLE}")

# create the list for f2py input
foreach ( file ${LIBSOURCES})
  set ( MODSOURCES ${MODSOURCES} ${PROJECT_SOURCE_DIR}/${file})
endforeach()

# find module name
execute_process( 
  # COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_config_var; print(get_config_var('EXT_SUFFIX'))"
  COMMAND ${Python3_EXECUTABLE} -c "from distutils.sysconfig import get_config_var; print(get_config_var('EXT_SUFFIX'))"
  OUTPUT_VARIABLE MODSUFFIX
  RESULT_VARIABLE MODSUFFIX_NOT_FOUND
  )
STRING(STRIP ${MODSUFFIX} MODSUFFIX)

if(MODNAME_NOT_FOUND)
  message(FATAL_ERROR "python3 EXT_SUFFIX not found")
endif()
message(STATUS "module fastpost${MODSUFFIX} will be build")

# add f2py target/command
add_custom_target(fastpost ALL DEPENDS fastpost${MODSUFFIX} )
if ( CMAKE_SYSTEM_NAME STREQUAL Linux OR 
     CMAKE_SYSTEM_NAME STREQUAL Darwin )
    add_custom_command(
      OUTPUT fastpost${MODSUFFIX}
      COMMAND  ${F2PY_EXECUTABLE} -m fastpost -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_BINARY_DIR} -L${PROJECT_BINARY_DIR} -lfastcore -c ${MODSOURCES}
      DEPENDS fastcore ${MODSOURCES}
    )
elseif(CMAKE_SYSTEM_NAME STREQUAL Windows)
    message( FATAL_ERROR "configuration for windows is still pending." )
endif()

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )
  # We are building as stand-alone project
  install (FILES ${PROJECT_BINARY_DIR}/fastpost${MODSUFFIX} DESTINATION ${PROJECT_BINARY_DIR}/../../pysimpp)

  # uninstall target
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

  # get the options for building tests (default no)
  option( DOTESTS "Build google.test" OFF )
  if ( DOTESTS )
      set( DOTESTS ${DOTESTS} CACHE BOOL "Build google.test" )
  endif()

  if ( DOTESTS )
    # add the test fwk
    set ( GTESTFWK ON )
    set ( BOOSTFWK OFF )
    set ( TESTS OFF )
    set ( EXAMPLES OFF )

    add_subdirectory( f90tw )

    # add test for fastpost core functionality
    F90TWTEST ( test_vector3d 
      CFILES tests/test_vector3d.cpp
      FPPFILES tests/test_vector3d.fpp
      OUTDIR "${CMAKE_CURRENT_BINARY_DIR}"
      EXE ${FPP_EXE} OPTIONS ${FPP_OPTIONS} SWITCH ${FPP_SWITCH}
      INCLUDES ${f90tw_INCLUDE_DIR}
      MODULES ${f90tw_MODULES_DIR} 
      LIBRARIES fastcore f90tw_gtest
      DEFINITIONS USEGTEST
    )

    F90TWTEST ( test_domain3d
      CFILES tests/test_domain3d.cpp
      F90FILES tests/test_domain3d.f90
      F2HFILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_domain3d.f90"
      OUTDIR "${CMAKE_CURRENT_BINARY_DIR}"
      EXE ${FPP_EXE} OPTIONS ${FPP_OPTIONS} SWITCH ${FPP_SWITCH}
      INCLUDES ${f90tw_INCLUDE_DIR}
      MODULES ${f90tw_MODULES_DIR} 
      LIBRARIES fastcore f90tw_gtest
      DEFINITIONS USEGTEST
    )

    enable_testing()
    add_test( test_vector3d test_vector3d )
    add_test( test_domain3d test_domain3d )

  endif()
  
else()
  # We are building as part of Main project
  install (FILES ${PROJECT_BINARY_DIR}/fastpost${MODSUFFIX} DESTINATION ${SKBUILD_PROJECT_NAME})

endif()
