cmake_minimum_required(VERSION 2.6)
project(ppi_networkit)

set(BENCHMARK_OMP ON CACHE BOOL "Use OpenMP parallelism in benchmark")

if (BENCHMARK_OMP)
    add_definitions(-DBENCHMARK_USE_OMP)
endif(BENCHMARK_OMP)

# set path for NetworKit
set(NETWORKIT_ROOT_DIR ${ppi_networkit_SOURCE_DIR}/ext/NetworKit CACHE PATH "The root dir of the NetworKit repo")
set(NETWORKIT_LIB_DIR ${NETWORKIT_ROOT_DIR} CACHE PATH "The library path for NetworKit")
set(NETWORKIT_SRC_DIR ${NETWORKIT_ROOT_DIR}/networkit/cpp CACHE PATH "The source
path for NetworKit headers")

# pull and compile NetworKit
add_custom_command(OUTPUT ${NETWORKIT_ROOT_DIR}/build.conf.example
                   COMMAND hg clone https://algohub.iti.kit.edu/parco/NetworKit/NetworKit ext/NetworKit
                   WORKING_DIRECTORY ${ppi_networkit_SOURCE_DIR}
                   COMMENT "Cloning NetworKit")

# TODO: make version a parameter
add_custom_command(OUTPUT ${NETWORKIT_ROOT_DIR}/build.conf
                   COMMAND hg checkout 3.3
                   COMMAND cp build.conf.example build.conf
                   DEPENDS ${NETWORKIT_ROOT_DIR}/build.conf.example
                   WORKING_DIRECTORY ${NETWORKIT_ROOT_DIR}
                   COMMENT "Checking out version 3.3")

add_custom_command(OUTPUT ${NETWORKIT_ROOT_DIR}/libNetworKit-Core-Opt.a
                   COMMAND scons --target=Core --optimize=Opt
                   DEPENDS ${NETWORKIT_ROOT_DIR}/build.conf
                   WORKING_DIRECTORY ${NETWORKIT_ROOT_DIR}
                   COMMENT "Building NetworKit core")

add_custom_target(networkit-lib
                  DEPENDS ${NETWORKIT_ROOT_DIR}/libNetworKit-Core-Opt.a)

# set include and link directories for NetworKit
link_directories(${NETWORKIT_ROOT_DIR})
include_directories(${NETWORKIT_ROOT_DIR}/networkit/cpp)

# set to use c++11 standard and openmp
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -Wall -O3 -msse2 -fopenmp -fPIC ${CMAKE_CXX_FLAGS}")

# Add the SQLiteC++ wrapper library
# Turn off annoying default options
set(SQLITECPP_BUILD_TESTS OFF)
set(SQLITECPP_BUILD_EXAMPLES OFF)
set(SQLITECPP_RUN_DOXYGEN OFF)
set(SQLITECPP_RUN_CPPCHECK OFF)
set(SQLITECPP_RUN_CPPLINT OFF)
set(SQLITECPP_INTERNAL_SQLITE OFF CACHE BOOL "use internal sqlite3" FORCE)
add_subdirectory(ext/SQLiteCpp)
include_directories("${PROJECT_SOURCE_DIR}/ext/SQLiteCpp/include")

# include my own files
include_directories(${ppi_networkit_SOURCE_DIR}/src)

# build ppi_networkit
add_subdirectory(src)

# benchmarks for tissue-specific betweenness
add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark ppi_networkit dl)

# tests executable
add_executable(tests tests.cpp)
target_link_libraries(tests ppi_networkit)
