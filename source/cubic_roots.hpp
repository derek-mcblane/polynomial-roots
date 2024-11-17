#pragma once

#include "quadratic_roots.hpp"
#include "small_integral_powers.hpp"

#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <complex>
#include <type_traits>

namespace dm::math {

namespace internal {

#if 0
template <typename RootPair>
struct CubicRoots
{
    using RealT = typename RootPair::RealT;
    using ComplexT = typename RootPair::ComplexT;

    RealT x1;
    RootPair pair;

    [[nodiscard]] constexpr std::array<ComplexT, 3> roots() const noexcept
    {
        return {ComplexT{x1, 0}, pair.r1(), pair.r2()};
    }

    [[nodiscard]] constexpr std::pair<std::array<RealT, 3>, std::size_t> real_roots() const noexcept
    {
        if constexpr (std::is_same_v<RootPair, RealRootPair<RootPair::RealT>>) {
            return {{x1, pair.x1(), pair.x2()}, 3};
        } else {
            return {{x1}, 1};
        }
    }
};
#endif

template <typename Real, template <typename> typename Complex = std::complex>
struct CubicRoots
{
    using RealT = Real;
    using ComplexT = Complex<Real>;

    Real x1;
    Real y1;
    Real x2;
    Real y2;
    Real x3;
    Real y3;

    std::array<ComplexT, 3> to_array()
    {
        return {ComplexT{x1, y1}, ComplexT{x2, y2}, ComplexT{x3, y3}};
    }
};

template <typename RealT>
struct CubicRealRoots
{
    using Real = RealT;

    Real x1;
    Real x2;
    Real x3;
    bool pair_real;

    std::pair<std::array<Real, 3>, std::size_t> to_array()
    {
        if (pair_real) {
            return {{x1, x2, x3}, 3};
        } else {
            return {{x1}, 1};
        }
    }
};

template <typename RealT>
struct MonicCubic
{
    using Real = RealT;

    /// x^3 + a[2]*x^2 + a[1]*x + a[0]
    std::array<Real, 3> a;

    [[nodiscard]] bool pair_real() const noexcept
    {
        return square(r()) <= -cube(q());
    }

    [[nodiscard]] CubicRoots<Real> roots() const noexcept
    {
        if (pair_real()) {
            return {three_x1(), 0, three_x2(), 0, three_x3(), 0};
        } else {
            return {one_x1(), 0, one_x2(), one_y2(), one_x2(), -one_y2()};
        }
    }

    [[nodiscard]] CubicRealRoots<Real> real_roots() const noexcept
    {
        CubicRealRoots<Real> roots;
        roots.pair_real = pair_real();
        if (roots.pair_real) {
            roots.x1 = three_x1();
            roots.x2 = three_x2();
            roots.x3 = three_x3();
        } else {
            roots.x1 = one_x1();
        }
        return roots;
    }

  private:
    [[nodiscard]] Real q() const noexcept
    {
        return a[1] / 3 - square(a[2]) / 9;
    }

    [[nodiscard]] Real r() const noexcept
    {
        return (a[1] * a[2] - 3 * a[0]) / 6 - cube(a[2]) / 27;
    }

    [[nodiscard]] Real A() const noexcept
    {
        return cbrt(abs(r()) + sqrt(square(r()) + cube(q())));
    }

    [[nodiscard]] Real one_t1() const noexcept
    {
        assert(!pair_real());
        return (r() >= 0) ? A() - q() / A() : q() / A() - A();
    }

    [[nodiscard]] Real one_x1() const noexcept
    {
        assert(!pair_real());
        return one_t1() - a[2] / 3;
    }

    [[nodiscard]] Real one_x2() const noexcept
    {
        assert(!pair_real());
        return -one_t1() / 2 - a[2] / 3;
    }

    [[nodiscard]] Real one_y2() const noexcept
    {
        assert(!pair_real());
        return sqrt(3) / 2 * (A() + q() / A());
    }

    [[nodiscard]] Real three_theta() const noexcept
    {
        assert(pair_real());
        return (q() != 0) ? acos(r() / sqrt(cube(-q()))) : 0;
    }

    [[nodiscard]] Real three_phi1() const noexcept
    {
        assert(pair_real());
        return three_theta() / 3;
    }

    [[nodiscard]] Real three_phi2() const noexcept
    {
        assert(pair_real());
        return three_phi1() - 2 * M_PI / 3;
    }

    [[nodiscard]] Real three_phi3() const noexcept
    {
        assert(pair_real());
        return three_phi1() + 2 * M_PI / 3;
    }

    [[nodiscard]] Real three_x1() const noexcept
    {
        assert(pair_real());
        return 2 * sqrt(-q()) * cos(three_phi1()) - a[2] / 3;
    }

    [[nodiscard]] Real three_x2() const noexcept
    {
        assert(pair_real());
        return 2 * sqrt(-q()) * cos(three_phi2()) - a[2] / 3;
    }

    [[nodiscard]] Real three_x3() const noexcept
    {
        assert(pair_real());
        return 2 * sqrt(-q()) * cos(three_phi3()) - a[2] / 3;
    }
};

} // namespace internal

template <typename Real, typename Coefficients>
[[nodiscard]] auto monic_cubic_roots(const Coefficients& c) noexcept -> std::array<std::complex<Real>, 3>
{
    return internal::MonicCubic<Real>{c[0], c[1], c[2]}.roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto cubic_roots(const Coefficients& c) noexcept -> std::array<std::complex<Real>, 3>
{
    return internal::MonicCubic<Real>{c[0] / c[3], c[1] / c[3], c[2] / c[3]}.roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto monic_cubic_real_roots(const Coefficients& c) noexcept -> std::pair<std::array<Real, 3>, std::size_t>
{
    return internal::MonicCubic<Real>{c[0], c[1], c[2]}.real_roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto cubic_real_roots(const Coefficients& c) noexcept -> std::pair<std::array<Real, 3>, std::size_t>
{
    if (c[3] == 0) {
        const auto [cubic_roots, n_roots] = quadratic_real_roots(std::array{c[0], c[1], c[2]});
        std::array<Real, 3> roots;
        std::copy_n(std::begin(cubic_roots), n_roots, std::begin(roots));
        return {roots, n_roots};
    }
    return internal::MonicCubic<Real>{c[0] / c[3], c[1] / c[3], c[2] / c[3]}.real_roots().to_array();
}

} // namespace dm::math
