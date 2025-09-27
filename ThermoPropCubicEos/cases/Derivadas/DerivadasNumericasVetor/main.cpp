#include "../../ThermoPropCubicEos.hpp"

double f1(double x){ return x*x;}
double fn1(double x, double dx){ return (f1(x+dx)- f1(x)/dx);}
double fa1(double x){ return 2;}

double f2(double x){ return x*x*x;}
double fn2(double x, double dx){ return (f2(x+dx)- f2(x)/dx);}
double fa2(double x){ return 3;}

double f3(double x){ return log(x);}
double fn3(double x, double dx){ return (f3(x+dx)- f3(x)/dx);}
double fa3(double x){ return 1/x;}

double f4(double x){ return exp(x);}
double fn4(double x, double dx){  return (f4(x+dx)- f4(x)/dx);}
double fa4(double x){ return exp(x);}

double f5(double x){ return cos(x);}
double fn5(double x, double dx){  return (f5(x+dx)- f5(x)/dx);}
double fa5(double x){ return -sin(x);}

int main ()
{
  auto x = 1.0;
  auto epsilon = std::numeric_limits<double>::epsilon();
  auto dx = epsilon;
  
  double fAnatical1 = fa1(x);
  double fAnatical2 = fa2(x);
  double fAnatical3 = fa3(x);
  double fAnatical4 = fa4(x);
  double fAnatical5 = fa5(x);

  //std::cout << std::fixed << std::setprecision(15);
  std::vector<double(*)(double)> f = {f1, f2, f3, f4, f5};
  std::vector<double(*)(double, double)> fn = {fn1, fn2, fn3, fn4, fn5};
  std::vector<double(*)(double)> fa  = {fa1, fa2, fa3, fa4, fa5};

  for (int i = 1; i <= 10; ++i) {
      double dx = std::pow(10.0, -i);
      for (size_t j = 0; j < f.size(); ++j) {
            double fNumerical = fn[j](x, dx);
            double fAnalitical = fa[j](x);
            //double erro = std::abs(fNumerical/ fAnalitical - 1.0);
            std::cout << "Erro da equação["<<j+1 <<"]com epsilon["<< dx <<"]: " << std::setprecision(15) << abs (fNumerical/fAnalitical - 1) << "\n" << std::endl;} 
  }
  return 0;
}


