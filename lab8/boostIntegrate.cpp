#include <fstream>
#include <vector>
#include <boost/numeric/odeint.hpp>


const double sigma = 10.0;
const double rho = 28.0;
const double beta = 8.0 / 3.0;

typedef std::vector<double> state;
std::ofstream outfile("boost.csv");

void lorenzSystem(const state &x, state &dxdt, double ){

  dxdt[0] = sigma * x[1] - sigma * x[0];
  dxdt[1] = -x[0] * x[2] + rho * x[0] - x[1];
  dxdt[2] = x[0] * x[1] - beta * x[2];
}

void saveCoord(const state &x, const double ){

  outfile << x[0] << "," << x[1] << "," << x[2] << std::endl;
}

int main(void){

  state x0 = {1.0, 1.0, 1.0};
  double t0 = 0.0, t = 25.0, dt = 0.1;
  boost::numeric::odeint::integrate(lorenzSystem, x0, t0, t, dt, saveCoord);
}

// compilation
// g++ -I 'C:/Program Files/boost_1_69_0' -o boostIntegrate boostIntegrate.cpp
// .\boostIntegrate