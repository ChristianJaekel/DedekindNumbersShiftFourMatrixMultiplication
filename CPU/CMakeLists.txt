cmake_minimum_required(VERSION 3.16)
project(dedekind)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-w -O3 -march=native -DEIGEN_NO_DEBUG")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# define DATA_DIR variable
set(DATA_DIR "${CMAKE_SOURCE_DIR}/../data/")

# configure header file with DATA_DIR value
configure_file(
  ${CMAKE_SOURCE_DIR}/src/DataDir.h.in
  ${CMAKE_SOURCE_DIR}/src/DataDir.h
)

# find Intel TBB
find_package(TBB REQUIRED)

# find Eigen
find_package(Eigen3 REQUIRED)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)

add_executable(
  dedekind
  ../src/FDL.cpp 
  ../src/IVLCalc.cpp         
  ../src/Timer.cpp 
  ../src/Uint128.cpp        
  src/main.cpp
  src/Enumerator.cpp
)

target_include_directories(dedekind PRIVATE src ../src ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(dedekind TBB::tbb Eigen3::Eigen)

add_subdirectory(benchmark)
add_subdirectory(test)
