#pragma once

#include <cstddef>
#include <utility>

namespace dm::math {

namespace internal {
template <typename Real, std::size_t... I>
Real ipow_(const Real value, std::index_sequence<I...>)
{
    return ((static_cast<void>(I), value) * ...);
}
} // namespace internal

template <std::size_t N, typename Real>
Real ipow(const Real value)
{
    return internal::ipow_(value, std::make_index_sequence<N>());
}

template <typename Real>
Real square(const Real value)
{
    return ipow<2>(value);
}

template <typename Real>
Real cube(const Real value)
{
    return ipow<3>(value);
}

} // namespace dm::math
