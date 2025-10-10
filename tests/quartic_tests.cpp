#include "quartic_roots.hpp"
#include "root_pair.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <type_traits>

using namespace dm::math;

using namespace testing;

template <typename RealType>
testing::internal::FloatingEqMatcher<RealType> RealNear(RealType rhs, RealType max_abs_error)
{
    return testing::internal::FloatingEqMatcher<RealType>(rhs, false, max_abs_error);
}

MATCHER_P2(
    ComplexNear,
    value,
    max_abs_error,
    "real " +
        DescribeMatcher<typename std::remove_reference_t<arg_type>::value_type>(
            RealNear(value.real(), max_abs_error), negation
        ) +
        ", imaginary " +
        DescribeMatcher<typename std::remove_reference_t<arg_type>::value_type>(
            RealNear(value.imag(), max_abs_error), negation
        )
)
{
    return ExplainMatchResult(RealNear(value.real(), max_abs_error), arg.real(), result_listener) &&
           ExplainMatchResult(RealNear(value.imag(), max_abs_error), arg.imag(), result_listener);
}

template <typename PairOne, typename PairTwo>
struct QuarticTestParams
{
    using PairOneT = PairOne;
    using PairTwoT = PairTwo;
    using RealT = typename PairOneT::RealT;
    using ComplexT = typename PairOneT::ComplexT;

    static_assert(
        std::is_same<typename PairOneT::RealT, typename PairTwoT::RealT>::value, "pair real types must match"
    );
    static_assert(
        std::is_same<typename PairOneT::ComplexT, typename PairTwoT::ComplexT>::value, "pair complex types must match"
    );

    constexpr QuarticTestParams(const PairOneT& p1, const PairTwoT& p2) : p1_(p1), p2_(p2) {}

    [[nodiscard]] RealT A3() const noexcept
    {
        return -(x1() + x2() + x3() + x4());
    }

    [[nodiscard]] RealT A2() const noexcept
    {
        return x1() * x2() + y1() * y1() + (x1() + x2()) * (x3() + x4()) + x3() * x4() + y3() * y3();
    }

    [[nodiscard]] RealT A1() const noexcept
    {
        return -((x1() * x2() + y1() * y1()) * (x3() + x4()) + (x3() * x4() + y3() * y3()) * (x1() + x2()));
    }

    [[nodiscard]] RealT A0() const noexcept
    {
        return (x1() * x2() + y1() * y1()) * (x3() * x4() + y3() * y3());
    }

    [[nodiscard]] ComplexT r1() const noexcept
    {
        return p1_.r1();
    }

    [[nodiscard]] ComplexT r2() const noexcept
    {
        return p1_.r2();
    }

    [[nodiscard]] ComplexT r3() const noexcept
    {
        return p2_.r1();
    }

    [[nodiscard]] ComplexT r4() const noexcept
    {
        return p2_.r2();
    }

  private:
    PairOne p1_;
    PairTwo p2_;

    [[nodiscard]] RealT x1() const noexcept
    {
        return p1_.x1();
    }

    [[nodiscard]] RealT y1() const noexcept
    {
        return p1_.y1();
    }

    [[nodiscard]] RealT x2() const noexcept
    {
        return p1_.x2();
    }

    [[nodiscard]] RealT y2() const noexcept
    {
        return p1_.y2();
    }

    [[nodiscard]] RealT x3() const noexcept
    {
        return p2_.x1();
    }

    [[nodiscard]] RealT y3() const noexcept
    {
        return p2_.y1();
    }

    [[nodiscard]] RealT x4() const noexcept
    {
        return p2_.x2();
    }

    [[nodiscard]] RealT y4() const noexcept
    {
        return p2_.y2();
    }
};

class QuarticPairOneRealTest
    : public ::testing::TestWithParam<QuarticTestParams<RealRootPair<double>, ComplexConjugateRootPair<double>>>
{};

static constexpr QuarticTestParams<RealRootPair<double>, ComplexConjugateRootPair<double>>
    quartic_pair_one_real_params[]{{RealRootPair<double>{10, 2}, ComplexConjugateRootPair<double>{3, 2}}};

INSTANTIATE_TEST_SUITE_P(Quartic, QuarticPairOneRealTest, testing::ValuesIn(quartic_pair_one_real_params));

TEST_P(QuarticPairOneRealTest, RootCheckEquations)
{
    const auto& p = GetParam();

    const auto roots = quartic_roots<double>(std::array{p.A0(), p.A1(), p.A2(), p.A3(), 1.0});

    const auto epsilon = 1e-9;
    EXPECT_THAT(
        roots,
        UnorderedElementsAre(
            ComplexNear(p.r1(), epsilon),
            ComplexNear(p.r2(), epsilon),
            ComplexNear(p.r3(), epsilon),
            ComplexNear(p.r4(), epsilon)
        )
    );
}

class QuarticPairTwoRealTest
    : public ::testing::TestWithParam<QuarticTestParams<ComplexConjugateRootPair<double>, RealRootPair<double>>>
{};

static constexpr QuarticTestParams<ComplexConjugateRootPair<double>, RealRootPair<double>>
    quartic_pair_two_real_params[]{{ComplexConjugateRootPair<double>{3, 2}, RealRootPair<double>{10, 2}}};

INSTANTIATE_TEST_SUITE_P(Quartic, QuarticPairTwoRealTest, testing::ValuesIn(quartic_pair_two_real_params));

TEST_P(QuarticPairTwoRealTest, RootCheckEquations)
{
    const auto& p = GetParam();

    const auto roots = quartic_roots<double>(std::array{p.A0(), p.A1(), p.A2(), p.A3(), 1.0});

    const auto epsilon = 1e-9;
    EXPECT_THAT(
        roots,
        UnorderedElementsAre(
            ComplexNear(p.r1(), epsilon),
            ComplexNear(p.r2(), epsilon),
            ComplexNear(p.r3(), epsilon),
            ComplexNear(p.r4(), epsilon)
        )
    );
}
