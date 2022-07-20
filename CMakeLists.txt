cmake_minimum_required(VERSION 3.16)
project(lbm)
project(lbm VERSION "0.0.2")

add_subdirectory(UnitTest)
add_subdirectory(perf)

configure_file(src/config.h.in config.h @ONLY)


set(CMAKE_CXX_STANDARD 17)

#todo: move to presets file
#include_directories(/home/svenb/build/omp411/include)
#link_directories(/home/svenb/build/omp411/lib)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

include_directories(src include)
include_directories(SYSTEM external)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(lbm
        src/main.cpp src/gridgenerator/gridGenerator.cpp src/lbm/solver.cpp
        src/gridgenerator/gridGenerator.h config.h src/globaltimers.h src/cartesiangrid.h
        src/gridcell_properties.h src/loadbalancing_weights.h src/geometry.h src/functions.h src/common/IO.h
        src/common/term.h src/interface/grid_interface.h src/gridgenerator/cartesiangrid_generation.h
        src/lbm/solver.h src/lbm/constants.h src/lbm/variables.h src/lbm/bnd/bnd.h src/common/surface.h
        src/lbm/solverExe.h src/common/configuration.h src/lbm/bnd/bnd_wall.h src/lbm/bnd/bnd_interface.h
        src/lbm/bnd/bnd_in_out.h src/lbm/bnd/bnd_periodic.h src/postprocessing.h src/postprocessing_cartesian.h
        src/lbm/lbm_postprocessing.h src/postprocessing_func.h src/common/line.h src/common/mesh.h
        src/lbm/bnd/bnd_dirichlet.h src/lbm/solver_inst_1d.cpp src/lbm/solver_inst_2d.cpp src/lbm/solver_inst_3d.cpp
        src/lbm/equilibrium_func.h
        src/lpt/particle.h src/lpt/forces.h src/lpt/solver.cpp src/lpt/solver.h src/common/sphere.h src/lpt/solverExe.h
        src/lpt/solver_inst_2d.cpp src/lpt/solver_inst_3d.cpp src/lpt/constants.h src/common/random.h
        src/lpt/analytical_solutions.h src/lbm/analytical_solutions.h src/lbm/bnd/bnd_neumann.h src/lbm/bnd/bnd_cell.h)

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
    add_custom_target(doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

