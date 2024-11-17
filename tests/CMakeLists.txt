include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE) # For Windows: Prevent overriding the parent project's compiler/linker settings
option(INSTALL_GMOCK "" OFF)
option(INSTALL_GTEST "" OFF)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_executable(CubicTests "")
target_sources(CubicTests PRIVATE cubic_tests.cpp)
target_include_directories(CubicTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(CubicTests PRIVATE PolynomialRoots gtest_main)
gtest_discover_tests(CubicTests)

add_executable(QuarticTests "")
target_sources(QuarticTests PRIVATE quartic_tests.cpp)
target_include_directories(QuarticTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(QuarticTests PRIVATE PolynomialRoots gtest_main)
gtest_discover_tests(QuarticTests)
