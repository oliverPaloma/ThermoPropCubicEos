#ifndef INCLUDES_H_INCLUDED
#define INCLUDES_H_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <assert.h>
#include <fstream>
#include <string.h>
#include <functional>

// The universal gas constant (in J/(mol*K))
const auto R = 8.3144621;

/// The value of NaN
constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();

using std::abs;
using std::log;
using std::sqrt;

using AlphaResult = std::tuple<double, double, double>;

/// Convenient alias for `std::function<R(Args...)>`.
template<typename F>
using Fn = std::function<F>;



#endif // INCLUDES_H_INCLUDED
