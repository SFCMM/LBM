project(UnitTest)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/src)

# openmpi
include_directories(/home/svenb/build/omp411/include)
link_directories(/home/svenb/build/omp411/lib)

# adding the Google_Tests_run target
add_executable(UnitTest test_hilbert.cpp test_math.cpp test_string_helper.cpp)
target_link_libraries(UnitTest gtest gtest_main gmock)

target_compile_options(UnitTest PUBLIC --std=c++17)