cmake_minimum_required(VERSION 3.25.0)

project(PolynomialRoots VERSION 0.1.0)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

set(PolynomialRoots_WARNING_OPTIONS
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wunused>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

add_compile_definitions(
    $<$<CONFIG:DEBUG>:POLYNOMIAL_ROOTS_DEBUG>
)

add_subdirectory(source)

option(PolynomialRoots_ENABLE_TESTING "Enable testing for PolynomialRoots" ON)
if (${PolynomialRoots_ENABLE_TESTING})
    enable_testing()
    add_subdirectory(tests)
endif()

export(EXPORT PolynomialRootsTargets
    NAMESPACE PolynomialRoots::
)
install(EXPORT PolynomialRootsTargets
    NAMESPACE PolynomialRoots::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PolynomialRoots
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/PolynomialRootsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion 
)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/PolynomialRootsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PolynomialRootsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PolynomialRoots
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/PolynomialRootsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PolynomialRootsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PolynomialRoots
)

export(PACKAGE PolynomialRoots)

include(CPack)
