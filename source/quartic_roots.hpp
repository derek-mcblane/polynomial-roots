#pragma once

#include "cubic_roots.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <utility>

namespace dm::math {

namespace internal {

template <typename Real, template <typename> typename Complex = std::complex>
struct QuarticRoots
{
    using RealT = Real;
    using ComplexT = Complex<Real>;

    Real x1;
    Real y1;
    Real x2;
    Real y2;
    Real x3;
    Real y3;
    Real x4;
    Real y4;

    std::array<ComplexT, 4> to_array()
    {
        return {ComplexT{x1, y1}, ComplexT{x2, y2}, ComplexT{x3, y3}, ComplexT{x4, y4}};
    }
};

template <typename Real>
struct QuarticRealRoots
{
    using RealT = Real;

    Real x1;
    Real x2;
    Real x3;
    Real x4;
    bool pair_one_real;
    bool pair_two_real;

    std::pair<std::array<Real, 4>, std::size_t> to_array()
    {
        if (pair_one_real && pair_two_real) {
            return {{x1, x2, x3, x4}, 4};
        } else if (pair_one_real) {
            return {{x1, x2}, 2};
        } else if (pair_two_real) {
            return {{x3, x4}, 2};
        } else {
            return {{}, 0};
        }
    }
};

template <typename Real, template <typename> typename Complex = std::complex>
class MonicQuartic
{
  public:
    using RealT = Real;
    using ComplexT = Complex<Real>;

    /// x^4 + A[3]*x^3 + A[2]*x^2 + A[1]*x + A[0]
    std::array<RealT, 4> A;

    [[nodiscard]] bool pair_one_real(const CubicRoots<RealT>& r) const noexcept
    {
        return radicand1(r) >= 0;
    }

    [[nodiscard]] bool pair_two_real(const CubicRoots<RealT>& r) const noexcept
    {
        return radicand2(r) >= 0;
    }

    [[nodiscard]] QuarticRoots<RealT> roots() const noexcept
    {
        QuarticRoots<RealT> roots;
        const auto r = resolvent_cubic_roots();
        if (pair_one_real(r)) {
            roots.x1 = sqrt(r.x1) - C() + sqrt(radicand1(r));
            roots.y1 = 0;
            roots.x2 = sqrt(r.x1) - C() - sqrt(radicand1(r));
            roots.y2 = 0;
        } else {
            roots.x1 = sqrt(r.x1) - C();
            roots.y1 = sqrt(-radicand1(r));
            roots.x2 = roots.x1;
            roots.y2 = -roots.y1;
        }
        if (pair_two_real(r)) {
            roots.x3 = -sqrt(r.x1) - C() + sqrt(radicand2(r));
            roots.y3 = 0;
            roots.x4 = -sqrt(r.x1) - C() - sqrt(radicand2(r));
            roots.y4 = 0;
        } else {
            roots.x3 = -sqrt(r.x1) - C();
            roots.y3 = sqrt(-radicand2(r));
            roots.x4 = roots.x3;
            roots.y4 = -roots.y3;
        }
        return roots;
    }

    [[nodiscard]] QuarticRealRoots<RealT> real_roots() const noexcept
    {
        const auto r = resolvent_cubic_roots();
        QuarticRealRoots<RealT> roots;
        roots.pair_one_real = pair_one_real(r);
        if (roots.pair_one_real) {
            roots.x1 = sqrt(r.x1) + sqrt(radicand1(r)) - C();
            roots.x2 = sqrt(r.x1) - sqrt(radicand1(r)) - C();
        }
        roots.pair_two_real = pair_two_real(r);
        if (roots.pair_two_real) {
            roots.x3 = -sqrt(r.x1) + sqrt(radicand2(r)) - C();
            roots.x4 = -sqrt(r.x1) - sqrt(radicand2(r)) - C();
        }
        return roots;
    }

  private:
    [[nodiscard]] RealT C() const noexcept
    {
        return A[3] / 4;
    }

    [[nodiscard]] RealT b0() const noexcept
    {
        return A[0] - A[1] * C() + A[2] * square(C()) - 3 * ipow<4>(C());
    }

    [[nodiscard]] RealT b1() const noexcept
    {
        return A[1] - 2 * A[2] * C() + 8 * cube(C());
    }

    [[nodiscard]] RealT b2() const noexcept
    {
        return A[2] - 6 * square(C());
    }

    [[nodiscard]] RealT sigma() const noexcept
    {
        return (b1() > 0) ? 1 : -1;
    }

    [[nodiscard]] CubicRoots<RealT> resolvent_cubic_roots() const noexcept
    {
        auto roots = MonicCubic<RealT>{-square(b1()) / 64, (square(b2()) - 4 * b0()) / 16, b2() / 2}.roots();
        if (roots.x1 < 0) {
            roots.x1 = 0;
        }
        if (roots.x2 * roots.x3 < 0) {
            if (roots.x2 > -roots.x3) {
                roots.x3 = 0;
            } else {
                roots.x2 = 0;
            }
        }
        return roots;
    }

    [[nodiscard]] RealT k(const CubicRoots<RealT>& r) const noexcept
    {
        return 2 * sigma() * sqrt(r.x2 * r.x3 + square(r.y2));
    }

    [[nodiscard]] RealT radicand1(const CubicRoots<RealT>& r) const noexcept
    {
        return r.x2 + r.x3 - k(r);
    }

    [[nodiscard]] RealT radicand2(const CubicRoots<RealT>& r) const noexcept
    {
        return r.x2 + r.x3 + k(r);
    }
};

} // namespace internal

template <typename Real, typename Coefficients>
[[nodiscard]] auto monic_quartic_roots(const Coefficients& c) -> std::array<std::complex<Real>, 4>
{
    return internal::MonicQuartic<Real>{c[0], c[1], c[2], c[3]}.roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto quartic_roots(const Coefficients& c) -> std::array<std::complex<Real>, 4>
{
    return internal::MonicQuartic<Real>{c[0] / c[4], c[1] / c[4], c[2] / c[4], c[3] / c[4]}.roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto monic_quartic_real_roots(const Coefficients& c) -> std::pair<std::array<Real, 4>, std::size_t>
{
    return internal::MonicQuartic<Real>{c[0], c[1], c[2], c[3]}.real_roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto quartic_real_roots(const Coefficients& c) -> std::pair<std::array<Real, 4>, std::size_t>
{
    if (c[4] == 0) {
        const auto [cubic_roots, n_roots] = cubic_real_roots<Real>(std::array{c[0], c[1], c[2], c[3]});
        std::array<Real, 4> roots;
        std::copy_n(std::begin(cubic_roots), n_roots, std::begin(roots));
        return {roots, n_roots};
    }
    return internal::MonicQuartic<Real>{c[0] / c[4], c[1] / c[4], c[2] / c[4], c[3] / c[4]}.real_roots().to_array();
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto monic_quartic_real_roots_sorted(const Coefficients& coefficients)
    -> std::pair<std::array<Real, 4>, std::size_t>
{
    auto [roots, n_real_roots] = monic_quartic_real_roots<Real>(coefficients);
    std::sort(roots.begin(), roots.begin() + n_real_roots);
    return {roots, n_real_roots};
}

template <typename Real, typename Coefficients>
[[nodiscard]] auto quartic_real_roots_sorted(const Coefficients& coefficients)
    -> std::pair<std::array<Real, 4>, std::size_t>
{
    auto [roots, n_real_roots] = quartic_real_roots<Real>(coefficients);
    std::sort(roots.begin(), roots.begin() + n_real_roots);
    return {roots, n_real_roots};
}

} // namespace dm::math
