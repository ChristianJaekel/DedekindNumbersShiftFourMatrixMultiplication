cmake_minimum_required(VERSION 3.16)
project(benchmark)

# download and include Google Benchmark
include(FetchContent)
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG main
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Suppressing benchmark's gtest tests" FORCE)
FetchContent_MakeAvailable(benchmark)

add_executable(
  dedekind_benchmarks
  main.cpp
  DedekindNumberBenchmark.cpp
  IntervalLengthBenchmark.cpp
  ../../src/FDL.cpp 
  ../../src/IVLCalc.cpp         
  ../../src/Timer.cpp 
  ../../src/Uint128.cpp  
  ../src/Enumerator.cpp
)

set(CMAKE_BUILD_TYPE Release)

target_include_directories(dedekind_benchmarks  PRIVATE ../src ../../src)

target_link_libraries(dedekind_benchmarks benchmark TBB::tbb Eigen3::Eigen)