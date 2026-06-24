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
#include <iomanip>
#include <functional>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <streambuf>
#include <Eigen/Core>
#include <unordered_map>

#include <numeric>
#include <sstream>
#include <cassert>


using namespace autodiff;

// machine epsilon square root
const auto SQRTepsilon = std::sqrt(std::numeric_limits<double>::epsilon());

/// The following definitions are based on the International System of Units(SI)
/// 9th edition 2019.

// The universal gas constant (in J/(mol*K))
const auto R = 8.31446261815324;
const auto Nav = 6.02214076e23;

/// The value of NaN
constexpr auto NaN = std::numeric_limits<real>::quiet_NaN();

using std::abs;
using std::log;
using std::sqrt;
using std::exp;

using AlphaResult = std::tuple<real, real, real>;

/// Convenient alias for `std::function<R(Args...)>`.
template<typename F>
using Fn = std::function<F>;

template<typename T>
using Vec = std::vector<T>;

using String = std::string;

#endif // INCLUDES_H_INCLUDED
