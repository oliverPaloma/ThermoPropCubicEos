#include "../../ThermoPropCubicEos.hpp"

double f(double x);
double f_numerica(double x, double dx);
double f_analitica (double x);

int main ()
{

  auto x = 1.0;
  //auto epsilon = 10^-3;
  auto epsilon = std::numeric_limits<double>::epsilon();
  auto dx = epsilon;
  //f(x);
  double fAnatical = f_analitica(x);

  for (int i = 1; i <= 10; ++i) {
    double dx = std::pow(10.0, -i);
    double fNum = f_numerica(x, dx);
    
    //std::cout << std::fixed << std::setprecision(15);
    std::cout << std ::setprecision(16)<< "Derivada numérica: " << fNum << std::endl;
    std::cout << "Erro entre Numérica e analitica: " << abs (fNum/fAnatical - 1) << "\n" << std::endl;
  }

  
  //auto fNum = f_numerica(x, dx)
  return 0;
}


double f(double x){
  return x*x + 2 * x + 1; 
}

double f_numerica(double x, double dx) {
    return (f(x + dx) - f(x)) / dx;
}

double f_analitica(double x){
  return 2*x + 2;
}





 