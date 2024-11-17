#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <utility>

namespace dm::math {

namespace internal {

template <typename Real, template <typename> typename Complex = std::complex>
struct QuadraticRoots
{
    using RealT = Real;
    using ComplexT = Complex<Real>;

    Real x1;
    Real x2;
    Real y1;

    std::array<ComplexT, 2> to_array()
    {
        return {ComplexT{x1, y1}, ComplexT{x2, -y1}};
    }
};

template <typename Real>
struct QuadraticRealRoots
{
    using RealT = Real;

    RealT x1;
    RealT x2;
    bool pair_real;

    std::pair<std::array<RealT, 2>, size_t> to_array()
    {
        if (pair_real) {
            return {{x1, x2}, 2};
        } else {
            return {{}, 0};
        }
    }
};

template <typename Real>
struct MonicQuadratic
{
    using RealT = Real;

    std::array<RealT, 2> c;

    [[nodiscard]] QuadraticRoots<RealT> roots() const noexcept
    {
        if (pair_real()) {
            return {two_x1(), 0, two_x2(), 0};
        } else {
            return {one_x1(), one_y1(), one_x1(), -one_y1()};
        }
    }

    [[nodiscard]] QuadraticRealRoots<RealT> real_roots() const noexcept
    {
        if (pair_real()) {
            return {two_x1(), two_x2(), true};
        } else {
            return {0, 0, false};
        }
    }

  private:
    [[nodiscard]] RealT two_x1() const noexcept
    {
        return (-c[1] + sqrt(D())) / 2;
    }

    [[nodiscard]] RealT two_x2() const noexcept
    {
        return (-c[1] - sqrt(D())) / 2;
    }

    [[nodiscard]] RealT one_x1() const noexcept
    {
        return -c[1] / 2;
    }

    [[nodiscard]] RealT one_y1() const noexcept
    {
        return sqrt(-D()) / 2;
    }

    [[nodiscard]] bool pair_real() const noexcept
    {
        return D() >= 0;
    }

    [[nodiscard]] RealT D() const noexcept
    {
        return c[1] * c[1] - 4 * c[0];
    }
};

} // namespace internal

template <typename Real, template <typename> typename Complex = std::complex>
std::array<Complex<Real>, 2> quadratic_roots(const std::array<Real, 3>& c)
{
    return internal::MonicQuadratic{c[0] / c[2], c[1] / c[2]}.roots().to_array();
}

template <typename Real>
std::pair<std::array<Real, 2>, size_t> quadratic_real_roots(const std::array<Real, 3>& c)
{
    if (c[2] == 0) {
        if (c[1] == 0) {
            return {{}, 0};
        }
        return {{0, -c[0] / c[1]}, 2};
    }
    return internal::MonicQuadratic<Real>{c[0] / c[2], c[1] / c[2]}.real_roots().to_array();
}

} // namespace dm::math
