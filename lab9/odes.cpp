#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

class IODE
	{
  public:
    virtual ~IODE() = default;

  public:
    virtual std::vector<double> getDerivatives(std::vector<double> values) = 0;

};


class LorenzSystem : public IODE
{
  public:
    LorenzSystem(double rho, double sigma, double beta) : 
      rho(rho), sigma(sigma), beta(beta) {}

  private: 
    std::vector<double> getDerivatives(std::vector<double> values){
      std::vector<double> result;
      result.push_back(sigma * values[1] - sigma * values[0]);
      result.push_back(-values[0] * values[2] + rho * values[0] - values[1]);
      result.push_back(values[0] * values[1] - beta * values[2]);
      return result;
    }

  private:
    double rho;
    double sigma;
    double beta;
};


class SimplePendulum : public IODE
{
  public:
    SimplePendulum(double l) : l(l) {}

  private: 
    std::vector<double> getDerivatives(std::vector<double> values){
      std::vector<double> result;
      result.push_back(-9.80665/l * sin(values[0]));
      return result;
    }

  private:
    double l;
};


class Spring : public IODE
{
  public:
    Spring(double m, double k) : 
    m(m), k(k) {}

  private: 
    std::vector<double> getDerivatives(std::vector<double> values){
      std::vector<double> result;
      result.push_back(-k * values[0] / m);
      return result;
    }

  private:
    double m;
    double k;
};


class LotkaVolterra : public IODE
{
  public:
    LotkaVolterra(double a, double b, double c, double d) : 
    a(a), b(b), c(c), d(d) {}

  private: 
    std::vector<double> getDerivatives(std::vector<double> values){
      std::vector<double> result;
      result.push_back(-b*values[0] + a*b*values[0]*values[1]);
      result.push_back(c*values[1] - d*values[0]*values[1]);
      return result;
    }

  private:
    double a;
    double b;
    double c;
    double d;
};


class Decay : public IODE
{
  public:
    Decay(double lambda) : lambda(lambda) {}

  private: 
    std::vector<double> getDerivatives(std::vector<double> values){
      std::vector<double> result;
      result.push_back(-lambda * values[0]);
      return result;
    }

  private:
    double lambda;
};


namespace ODESolver
{

	class Lab8
  {
    public:
      Lab8(IODE& ode, double dt, int max_steps) : 
        ode(ode), dt(dt), max_steps(max_steps)  {}

      void eulerMethod(std::vector<double> values){

        int step = 1;
        std::vector<double> result;
        std::ofstream outfile("explicit_euler.csv");

        while(step <= max_steps){

          result = ode.getDerivatives(values);

          for(int i=0; i<result.size(); i++){
            values[i] += result[i] * dt;
          }
          
          for(int i=0; i<values.size(); i++){
            outfile << values[i];
            if(i != values.size()-1) outfile << ",";
            else outfile << std::endl;
          }

          step++;
        }

      }

      void backwardEulerMethod(std::vector<double> values){

        double err = 1e-6;
        int step = 1;
        std::vector<double> result, values0(values.size()), values1(values.size());
        
        std::ofstream outfile("implicit_euler.csv");

        while(step <= max_steps){

          result = ode.getDerivatives(values);
          for(int i=0; i<result.size(); i++){
            values0[i] = values0[i] + result[i] * dt;
          }

          result = ode.getDerivatives(values0);
          for(int i=0; i<result.size(); i++){
            values1[i] = values1[i] + result[i] * dt;
          }

          while(!isClose(values0, values1, err)){

            values0 = values1;
            
            result = ode.getDerivatives(values0);
            for(int i=0; i<result.size(); i++){
              values1[i] = values[0] + result[i] * dt;
            }
          }

          values = values1;

          for(int i=0; i<values.size(); i++){
            outfile << values[i];
            if(i != values.size()-1) outfile << ",";
            else outfile << std::endl;
          }

          step++;
        }
        
      }


      void midpointMethod(std::vector<double> values){

        int step = 1;
        std::vector<double> result, tmp(values.size());

        std::ofstream outfile("midpoint.csv");

        while(step <= max_steps){

          result = ode.getDerivatives(values);
          for(int i=0; i<result.size(); i++){
            tmp[i] = values[i] + result[i] * dt/2;
          }

          result = ode.getDerivatives(values);
          for(int i=0; i<result.size(); i++){
            values[i] = values[i] + result[i] * dt;
          }

          for(int i=0; i<values.size(); i++){
            outfile << values[i];
            if(i != values.size()-1) outfile << ",";
            else outfile << std::endl;
          }

          step++;
        }
        
      }


      void rungeKuttaMethod(std::vector<double> values){

        int step = 1;
        std::vector<double> k1, k2, k3, k4, tmp(values.size());
        
        std::ofstream outfile("runge_kutta.csv");

        while(step <= max_steps){

          k1 = ode.getDerivatives(values);
          for(int i=0; i<k1.size(); i++){
            k1[i] *= dt;
            tmp[i] = values[i] + k1[i]/2;
          }

          k2 = ode.getDerivatives(tmp);
          for(int i=0; i<k2.size(); i++){
            k2[i] *= dt;
            tmp[i] = values[i] + k2[i]/2;
          }

          k3 = ode.getDerivatives(tmp);
          for(int i=0; i<k3.size(); i++){
            k3[i] *= dt;
            tmp[i] = values[i] + k3[i];
          }

          k4 = ode.getDerivatives(tmp);
          for(int i=0; i<k4.size(); i++){
            k4[i] *= dt;
          }

          for(int i=0; i<values.size(); i++){
            values[i] = values[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
          }

          for(int i=0; i<values.size(); i++){
            outfile << values[i];
            if(i != values.size()-1) outfile << ",";
            else outfile << std::endl;
          }

          step++;
        }
        
      }
      

    private: 
      bool isClose(std::vector<double> v1, std::vector<double> v2, double err){
        for(int i=0; i<v1.size(); i++){
          if(abs(v1[i]-v2[i]) >= err) return false;
        }
        return true;
      }

    private:
      IODE& ode;
      double dt;
      int max_steps;
  };

  class Leapfrog
  {
    
    public:
      Leapfrog(IODE& ode, double dt, int max_steps) : 
          ode(ode), dt(dt), max_steps(max_steps)  {}

      void calculate(std::vector<double> x, std::vector<double> v){
        
        int step = 1;
        std::vector<double> a = ode.getDerivatives(x);
        std::vector<double> tmp(a.size());
        std::ofstream outfile("leapfrog.csv");

        while(step <= max_steps){

          for(int i=0; i<x.size(); i++){
            x[i] = x[i] + v[i] * dt + (a[i] * dt * dt)/2;
          }

          tmp = a;
          a = ode.getDerivatives(x);

          for(int i=0; i<v.size(); i++){
            v[i] = v[i] + (tmp[i] + a[i])*dt/2;
          }

          // format: t,x1,v1,a1,x2,v2,a2...
          outfile << dt*step << ",";
          for(int i=0; i<x.size(); i++){
              outfile << x[i] << "," << v[i] << "," << a[i];
              if(i != x.size()-1) outfile << ",";
              else outfile << std::endl;
            }
          
          step++;
        }

    }


    private:
      IODE& ode;
      double dt;
      int max_steps;
  };

}

void testSimplePendulum(){

  std::vector<double> initialX = {0.069813}; // 4 degrees
  //std::vector<double> initialX = {0.785398}; // 45 degrees
  std::vector<double> initialV = {0.0};
  SimplePendulum pendulum = SimplePendulum(1);

  std::unique_ptr<ODESolver::Leapfrog> solver = 
					std::make_unique<ODESolver::Leapfrog>(ODESolver::Leapfrog(pendulum, 0.01, 1000));

  solver->calculate(initialX, initialV);

}

void testSpring(){

  std::vector<double> initialX = {0.6};
  std::vector<double> initialV = {0.0};
  Spring spring = Spring(2, 300);

  std::unique_ptr<ODESolver::Leapfrog> solver = 
					std::make_unique<ODESolver::Leapfrog>(ODESolver::Leapfrog(spring, 0.01, 500));

  solver->calculate(initialX, initialV);

}

void testLV(){

  std::vector<double> initialValues = {10.0, 10.0};
  LotkaVolterra lv = LotkaVolterra(0.3, 0.4, 0.9, 0.4);

  std::unique_ptr<ODESolver::Lab8> solver = 
					std::make_unique<ODESolver::Lab8>(ODESolver::Lab8(lv, 0.01, 6000));

  solver->rungeKuttaMethod(initialValues);

}

void testDecay(){

  std::vector<double> initialValue = {1.0};
  Decay decay = Decay(1.0);

  std::unique_ptr<ODESolver::Lab8> solver = 
					std::make_unique<ODESolver::Lab8>(ODESolver::Lab8(decay, 0.01, 2500));

  solver->rungeKuttaMethod(initialValue);

}


int main(void){

  //testSimplePendulum();

  //testSpring();

  //testLV();

  testDecay();

  return 0;
}