#pragma once
#include "../includes.hpp"
#include "correlations.hpp"
#include "Lumping.hpp"
#include "cubic_eos_types.hpp"



// C++ includes
#include <array>
#include <complex>
#include <functional>
#include <assert.h>

/// Define a type that describes the roots of a cubic equation
template<typename T>
using CubicRoots = std::array<std::complex<T>, 3>;

/// Calculate the roots of a cubic equation using Cardano's method.
/// The calculation uses the approach presented in:
///
///    Nickalls, R. W. D. (2012). A New Approach to Solving the Cubic: Cardan’s
///    Solution Revealed. The Mathematical Gazette, 77(480), 354–359* to solve
///    the cubic equation: \f$ x^{3}+bx^{2}+cx+d=0 \f$.
///
/// @param b The coefficient *b* of the cubic equation
/// @param c The coefficient *c* of the cubic equation
/// @param d The coefficient *d* of the cubic equation
/// @return The three roots \f$ (r_1, r_2, r_3) \f$ that solve the cubic
/// equation. No ordering is performed on the roots!
template<typename T>
auto cardano(const T& b, const T& c, const T& d) -> CubicRoots<T>
{
    using std::abs;
    using std::acos;
    using std::cbrt;
    using std::cos;
    using std::sqrt;
    using std::complex;

    const auto xn = -b/3;
    const auto yn = d + xn*(c + xn*(b + xn));

    const auto delta2 = (b*b - 3*c)/9;

    const auto u = yn*yn;
    const auto v = 4*delta2*delta2*delta2;

    const auto discr = u - v;

    const auto eps = 100*std::numeric_limits<T>::epsilon();

    if(discr > eps) {
        const auto sqrtdiscr = sqrt(discr);
        const auto aux1 = cbrt( 0.5*(-yn + sqrtdiscr) );
        const auto aux2 = cbrt( 0.5*(-yn - sqrtdiscr) );
        const auto alpha = xn + aux1 + aux2;
        const auto z1 = alpha - xn;
        const auto aux3 = 0.5*sqrt(3*(z1*z1 - 4*delta2));
        const auto aux4 = xn - 0.5*z1;
        const auto beta = complex{aux4, -aux3};
        const auto gamma = complex{aux4, aux3};
        return {alpha, beta, gamma};
    }
    else if(discr < -eps) {
        const auto pi = 3.14159265358979323846;
        const auto delta = sqrt(delta2);
        const auto h = 2*delta2*delta;
        const auto theta = acos(-yn/h) / 3.0;
        const auto alpha = xn + 2*delta*cos(theta);
        const auto beta  = xn + 2*delta*cos(2*pi/3 - theta);
        const auto gamma = xn + 2*delta*cos(2*pi/3 + theta);
        return {gamma, beta, alpha};
    }
    else {
        const auto delta = cbrt( 0.5*yn );
        const auto alpha = xn + delta;
        const auto gamma = xn - 2*delta;
        if(alpha < gamma) return {alpha, alpha, gamma};
        else return {gamma, alpha, alpha};
    }
}

/// Calculate the root of a non-linear function using Newton's method.
/// @param f The function that returns a pair of \f$ f(x) \f$ and \f$ f^{\prime}(x) \f$.
/// @param x0 The initial guess for the iterative root calculation.
/// @param epsilon The tolerance used in \f$ |f(x)| < \epsilon \f$ to check convergence.
/// @param maxiter The maximum number of iterations.
/// @return A tuple containing the root, the number of iterations, and a boolean flag indicating success if true.
template<typename T>
auto newton(const std::function<std::tuple<T, T>(const T&)>& f, const T& x0, const T& epsilon, std::size_t maxiter) -> std::tuple<T, std::size_t, bool>
{
    using std::abs;
    assert(epsilon > 0.0);
    assert(maxiter > 0);
    T x = x0;
    for(auto i = 0; i < maxiter; ++i)
    {
        const auto [fx, dfx] = f(x);
        assert(dfx != 0.0);
        x -= fx/dfx;
        if(abs(fx) < epsilon)
            return { x, i, true };
    }
    return { x, maxiter, false };
}

/// Return all real roots of a group of roots
/// @param roots CubicRoots with of complex and real roots
/// @return A vector with all real roots
template<typename T>
auto realRoots(const CubicRoots<T>& roots) -> std::vector<T>
{
    std::vector<T> real_roots;
    if(std::get<0>(roots).imag() == 0.0)
        real_roots.push_back(std::get<0>(roots).real());
    if(std::get<1>(roots).imag() == 0.0)
        real_roots.push_back(std::get<1>(roots).real());
    if(std::get<2>(roots).imag() == 0.0)
        real_roots.push_back(std::get<2>(roots).real());

    return real_roots;
}

// Calculates the norm of a vector
template<typename T>
auto norm(Vec<T> const& v, int flag) -> T
{
  if (v.empty()) return 0.;
  static real sum = 0.;
  auto n = v.size();
  
  switch (flag)
  {
    case 0:  // Maximum norm
      static Vec<T> vaux;
      vaux.resize(n);
      for (auto i = 0; i < n; i++) vaux[i] = abs(v[i]);
      return *std::max_element(vaux.begin(),vaux.end());
      break;
    case 1:  // Taxicab Norm or Manhattan norm
      for (auto i = 0; i < n; i++) sum += abs(v[i]);
      return sum;
      break;
    default: // Euclidian Norm
      for (auto i = 0; i < n; i++) sum += v[i]*v[i];
      return sqrt(sum);
  }
}  

template<typename T>
auto lubksb(Vec<T> &a, int n, Vec<int> &idx, Vec<T> &b) -> void
{
	int i,ii=0,ip,j;
	T sum;

	for (i=0;i<n;i++) {
		ip=idx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
      for (j=ii-1;j<i;j++) sum -= a[j+n*i]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
 		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[j+n*i]*b[j];
		b[i]=sum/a[i+n*i];
	}
  
}

template<typename T>
auto ludcmp(Vec<T> &a, int n, Vec<int> &idx, T &d) -> void
{
	const T TINY = 1.0e-20;
	int i,imax,j,k;
	T big,dum,sum,temp;

	static Vec<T> vv; // different values of n can be used in the same run
  vv.resize(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=abs(a[j+n*i])) > big) big=temp;
		if (big == 0.0) std::cout << "Singular matrix in routine ludcmp" << std::endl; 
    vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[j+n*i];
			for (k=0;k<i;k++) sum -= a[k+n*i]*a[j+n*k];
			a[j+n*i]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[j+n*i];
			for (k=0;k<j;k++) sum -= a[k+n*i]*a[j+n*k];
			a[j+n*i]=sum;
			if ((dum=vv[i]*abs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[k+n*imax];
				a[k+n*imax]=a[k+n*j];
				a[k+n*j]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		idx[j]=imax;
		if (a[j+n*j] == 0.0) a[j+n*j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j+n*j]);
			for (i=j+1;i<n;i++) a[j+n*i] *= dum;
		}
	}
  
}

// Calculates mixed error
template<typename T>
auto mixed_error(Vec<T> const& deltaX, Vec<T> const& X, T const& absE) -> T
{
  auto n = X.size();
  assert(("Calculating mixed error: deltaX and X should be non-empty vectors",(!deltaX.empty() && !X.empty())));
  assert(("Calculating mixed error: deltaX and X should have the same size",(deltaX.size() == n)));

  static real sqr_sum;
  sqr_sum = 0.;
  static real relE = absE*100.;
//  for (auto i = 0; i < n; i++) sqr_sum += (std::abs(deltaX[i]) / (std::abs(X[i])*relE + absE))* 
//                                     (std::abs(deltaX[i]) / (std::abs(X[i])*relE + absE));
  for (auto i = 0; i < n; i++) sqr_sum += (deltaX[i] / (X[i]*relE + absE))* 
                                      (deltaX[i] / (X[i]*relE + absE));
  return sqrt(sqr_sum/n);
} 

// Function that returns the index of the maximum absolute value of a vector
template<typename T>
auto idxMaxAbsVec(Vec<T> &v) -> int
{
  if (v.empty()) return 0.;
  for (auto i = 0; i < v.size(); i++) v[i] = abs(v[i]);
  return std::max_element(v.begin(),v.end()) - v.begin();
}  