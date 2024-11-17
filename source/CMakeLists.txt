add_library(PolynomialRoots INTERFACE)
add_library(PolynomialRoots::PolynomialRoots ALIAS PolynomialRoots)
target_compile_features(PolynomialRoots INTERFACE cxx_std_11)
target_compile_options(PolynomialRoots INTERFACE ${PolynomialRoots_WARNING_OPTIONS})
target_include_directories(PolynomialRoots INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_sources(PolynomialRoots
PUBLIC FILE_SET HEADERS FILES
    small_integral_powers.hpp
    root_pair.hpp
    quadratic_roots.hpp
    cubic_roots.hpp
    quartic_roots.hpp
)
install(TARGETS PolynomialRoots EXPORT PolynomialRootsTargets
    FILE_SET HEADERS
    INCLUDES DESTINATIION ${CMAKE_INSTALL_INCLUDEDIR}
)
