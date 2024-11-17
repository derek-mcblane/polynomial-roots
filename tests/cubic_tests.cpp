#include <gtest/gtest.h>

#include "cubic_roots.hpp"

template <typename Real>
struct ThreeRealRootCubicTestParams
{
    Real x1;
    Real x2;
    Real x3;

    Real a2() const
    {
        return -(x1 + x2 + x3);
    }

    Real a1() const
    {
        return x1 * (x2 + x3) + x2 * x3;
    }

    Real a0() const
    {
        return -x1 * x2 * x3;
    }
};

class ThreeRealRootCubicTest : public ::testing::TestWithParam<ThreeRealRootCubicTestParams<double>>
{};

static constexpr ThreeRealRootCubicTestParams<double> three_real_cubic_params[]{
    {3, 2, 1}, {2, 1, -10}, {1000.0, 0.5, -100.0}
};
INSTANTIATE_TEST_SUITE_P(Cubic, ThreeRealRootCubicTest, testing::ValuesIn(three_real_cubic_params));

TEST_P(ThreeRealRootCubicTest, RootCheckEquations)
{
    const auto& p = GetParam();
    ASSERT_GE(p.x1, p.x2);
    ASSERT_GE(p.x2, p.x3);

    const auto roots = dm::math::cubic_roots<double>(std::array{p.a0(), p.a1(), p.a2(), 1.0});

    const auto epsilon = 1e-9;
    EXPECT_NEAR(roots[0].real(), p.x1, epsilon);
    EXPECT_NEAR(roots[1].real(), p.x2, epsilon);
    EXPECT_NEAR(roots[2].real(), p.x3, epsilon);
    EXPECT_NEAR(roots[0].imag(), 0.0, epsilon);
    EXPECT_NEAR(roots[1].imag(), 0.0, epsilon);
    EXPECT_NEAR(roots[2].imag(), 0.0, epsilon);
}

TEST_P(ThreeRealRootCubicTest, RealRootCheckEquations)
{
    const auto& p = GetParam();
    ASSERT_GE(p.x1, p.x2);
    ASSERT_GE(p.x2, p.x3);

    const auto [roots, n_roots] = dm::math::cubic_real_roots<double>(std::array{p.a0(), p.a1(), p.a2(), 1.0});

    EXPECT_EQ(n_roots, 3);

    const auto epsilon = 1e-9;
    EXPECT_NEAR(roots[0], p.x1, epsilon);
    EXPECT_NEAR(roots[1], p.x2, epsilon);
    EXPECT_NEAR(roots[2], p.x3, epsilon);
}

template <typename Real>
struct OneRealRootCubicTestParams
{
    Real x1;
    Real x2;
    Real y2;

    Real a2() const
    {
        return -(x1 + x2 + x2);
    }
    Real a1() const
    {
        return 2 * x1 * x2 + x2 * x2 + y2 * y2;
    }
    Real a0() const
    {
        return -x1 * (x2 * x2 + y2 * y2);
    }
};

class OneRealRootCubicTest : public ::testing::TestWithParam<OneRealRootCubicTestParams<double>>
{};

static constexpr OneRealRootCubicTestParams<double> one_real_cubic_params[]{
    {3, 2, 1}, {1, -10, 2}, {0.5, -100.0, 1000.0}
};
INSTANTIATE_TEST_SUITE_P(Cubic, OneRealRootCubicTest, testing::ValuesIn(one_real_cubic_params));

TEST_P(OneRealRootCubicTest, RootCheckEquations)
{
    const auto& p = GetParam();

    if (p.y2 != 0) {
        ASSERT_NE(p.a1(), 2 * p.x1 * p.x2 + p.x2 * p.x2) << "y2 is too small relative to x1 and x2";
    }

    const auto roots = dm::math::cubic_roots<double>(std::array{p.a0(), p.a1(), p.a2(), 1.0});

    const auto epsilon = 1e-9;
    EXPECT_NEAR(roots[0].real(), p.x1, epsilon);
    EXPECT_NEAR(roots[1].real(), p.x2, epsilon);
    EXPECT_NEAR(roots[2].real(), p.x2, epsilon);
    EXPECT_NEAR(roots[0].imag(), 0.0, epsilon);
    if (p.y2 > 0) {
        EXPECT_NEAR(roots[1].imag(), p.y2, epsilon);
        EXPECT_NEAR(roots[2].imag(), -p.y2, epsilon);
    } else {
        EXPECT_NEAR(roots[1].imag(), -p.y2, epsilon);
        EXPECT_NEAR(roots[2].imag(), p.y2, epsilon);
    }
}

TEST_P(OneRealRootCubicTest, RealRootCheckEquations)
{
    const auto& p = GetParam();

    if (p.y2 != 0) {
        ASSERT_NE(p.a1(), 2 * p.x1 * p.x2 + p.x2 * p.x2) << "y2 is too small relative to x1 and x2";
    }

    const auto [roots, n_roots] = dm::math::cubic_real_roots<double>(std::array{p.a0(), p.a1(), p.a2(), 1.0});

    EXPECT_EQ(n_roots, 1);

    const auto epsilon = 1e-9;
    EXPECT_NEAR(roots[0], p.x1, epsilon);
}
