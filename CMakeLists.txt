cmake_minimum_required(VERSION 3.16)
project(lbm)
project(lbm VERSION "0.0.1")

add_subdirectory(UnitTest)
add_subdirectory(perf)

configure_file(src/config.h.in config.h @ONLY)


set(CMAKE_CXX_STANDARD 17)

#todo: move to presets file
include_directories(/home/svenb/build/omp411/include)
link_directories(/home/svenb/build/omp411/lib)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

include_directories(external src include)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(lbm
        src/main.cpp src/gridgenerator/gridGenerator.cpp src/lbm_solver.cpp
        src/gridgenerator/gridGenerator.h config.h
        src/globaltimers.h src/cartesiangrid.h
        src/gridcell_properties.h src/WeightMethod.h src/geometry.h src/functions.h src/common/IO.h src/term.h
        src/interface/grid_interface.h src/gridgenerator/cartesiangrid_generation.h src/lbm_solver.h
        src/lbm_constants.h src/pv.h src/interface/lbm_interface.h src/lbm_bnd.h src/common/surface.h src/lbm_solverExe.h src/configuration.h src/lbm_bnd_wall.h src/lbm_bnd_interface.h src/lbm_bnd_in_out.h)

#TARGET_LINK_LIBRARIES(lbm PUBLIC mpi)
target_link_libraries(lbm PUBLIC MPI::MPI_CXX dl tinfo pthread)

option(BUILD_DOC "Build documentation" OFF)
# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

