#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <algorithm>

namespace LorenzSystem
{

	class Solver
  {
    public:
      Solver(double rho, double sigma, double beta, double dt, int max_steps) : 
        rho(rho), sigma(sigma), beta(beta), dt(dt), max_steps(max_steps)  {}

      void eulerMethod(double x, double y, double z){

        int step = 1;
        std::vector<double> result;
        std::ofstream outfile("explicit_euler.csv");

        while(step <= max_steps){

          result = getValues(x, y, z);
          x = x + result[0] * dt;
          y = y + result[1] * dt;
          z = z + result[2] * dt;
          
          outfile << x << "," << y << "," << z << std::endl;

          step++;
        }

      }

      void backwardEulerMethod(double x, double y, double z){

        double err = 1e-6;
        int step = 1;
        std::vector<double> result;
        std::ofstream outfile("implicit_euler.csv");

        while(step <= max_steps){

          result = getValues(x, y, z);
          double x0 = x + dt * result[0];
          double y0 = y + dt * result[1];
          double z0 = z + dt * result[2];

          result = getValues(x0, y0, z0);
          double x1 = x + dt * result[0];
          double y1 = y + dt * result[1];
          double z1 = z + dt * result[2]; 

          while(abs(x0-x1) >= err || abs(y0-y1) >= err || abs(z0-z1) >= err){

            x0 = x1;
            y0 = y1;
            z0 = z1;

            result = getValues(x0, y0, z0);
            double x1 = x + dt * result[0];
            double y1 = y + dt * result[1];
            double z1 = z + dt * result[2];
          }

          x = x1;
          y = y1;
          z = z1;

          outfile << x << "," << y << "," << z << std::endl;

          step++;
        }
        
      }


      void midpointMethod(double x, double y, double z){

        int step = 1;
        std::vector<double> result;
        std::ofstream outfile("midpoint.csv");

        while(step <= max_steps){

          result = getValues(x, y, z);
          double tmpx = x + dt/2 * result[0];
          double tmpy = y + dt/2 * result[1];
          double tmpz = z + dt/2 * result[2];

          result = getValues(tmpx, tmpy, tmpz);
          x = x + result[0] * dt;
          y = y + result[1] * dt;
          z = z + result[2] * dt;

          outfile << x << "," << y << "," << z << std::endl;

          step++;
        }
        
      }


      void rungeKuttaMethod(double x, double y, double z){

        std::vector<double> k1, k2, k3, k4;
        int step = 1;
        std::ofstream outfile("runge_kutta.csv");

        while(step <= max_steps){

          k1 = getValues(x, y, z);
          for(int i=0; i<3; i++) k1[i] *= dt;

          k2 = getValues(x+k1[0]/2, y+k1[1]/2, z+k1[2]/2);
          for(int i=0; i<3; i++) k2[i] *= dt;

          k3 = getValues(x+k2[0]/2, y+k2[1]/2, z+k2[2]/2);
          for(int i=0; i<3; i++) k3[i] *= dt;

          k4 = getValues(x+k3[0], y+k3[1], z+k3[2]);
          for(int i=0; i<3; i++) k4[i] *= dt;

          x = x + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6;
          y = y + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6;
          z = z + (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6;


          outfile << x << "," << y << "," << z << std::endl;

          step++;
        }
        
      }
      

    private: 
      std::vector<double> getValues(double x, double y, double z){
        std::vector<double> result;
        result.push_back(sigma * y - sigma * x);
        result.push_back(-x * z + rho * x - y);
        result.push_back(x * y - beta * z);
        return result;
      }

    private:
      double rho;
      double sigma;
      double beta;
      double dt;
      int max_steps;
  };

}


int main(void){

  std::unique_ptr<LorenzSystem::Solver> lorenz = 
					std::make_unique<LorenzSystem::Solver>(LorenzSystem::Solver(28, 10, 8/3, 0.01, 2500));

  lorenz->eulerMethod(1, 1, 1);

  lorenz->backwardEulerMethod(1, 1, 1);

  lorenz->midpointMethod(1, 1, 1);

  lorenz->rungeKuttaMethod(1, 1, 1);

  return 0;
}