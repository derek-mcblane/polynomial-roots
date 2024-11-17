#pragma once

#include <array>
#include <complex>
#include <type_traits>
#include <utility>

namespace dm::math {

template <typename Real, template <typename> typename Complex = std::complex>
class ComplexConjugateRootPair
{
  public:
    using RealT = Real;
    using ComplexT = Complex<Real>;

    constexpr ComplexConjugateRootPair(RealT x1, RealT y1) : x1_(x1), y1_(y1) {}

    RealT x1() const
    {
        return x1_;
    }

    RealT x2() const
    {
        return x1_;
    }

    RealT y1() const
    {
        return y1_;
    }

    RealT y2() const
    {
        return -y1_;
    }

    ComplexT r1() const
    {
        return {x1(), y1()};
    }

    ComplexT r2() const
    {
        return {x2(), y2()};
    }

  private:
    RealT x1_;
    RealT y1_;
};

template <typename Real, template <typename> typename Complex = std::complex>
class RealRootPair
{
  public:
    using RealT = Real;
    using ComplexT = Complex<Real>;

    constexpr RealRootPair(RealT x1, RealT x2) noexcept : x1_(x1), x2_(x2) {}

    [[nodiscard]] constexpr RealT x1() const noexcept
    {
        return x1_;
    }

    [[nodiscard]] constexpr RealT x2() const noexcept
    {
        return x2_;
    }

    [[nodiscard]] constexpr RealT y1() const noexcept
    {
        return 0;
    }

    [[nodiscard]] constexpr RealT y2() const noexcept
    {
        return 0;
    }

    [[nodiscard]] constexpr ComplexT r1() const noexcept
    {
        return {x1(), y1()};
    }

    [[nodiscard]] constexpr ComplexT r2() const noexcept
    {
        return {x2(), y2()};
    }

  private:
    RealT x1_;
    RealT x2_;
};

} // namespace dm::math
