# download and include Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.11.0
)
FetchContent_MakeAvailable(googletest)

add_executable(
    dedekind_tests
    main.cpp
    DedekindNumberTest.cpp
    IntervalLengthTest.cpp
    ../../src/FDL.cpp 
    ../../src/IVLCalc.cpp         
    ../../src/Timer.cpp 
    ../../src/Uint128.cpp  
    ../src/Enumerator.cpp
)

target_include_directories(dedekind_tests  PRIVATE ../src ../../src)

target_link_libraries(dedekind_tests PUBLIC gtest_main Eigen3::Eigen TBB::tbb)