#include <iostream>
#include <functional>
#include <memory>
#include <math.h>
#include <time.h>

namespace integration
{

	class IIntegration
	{
    public:
      virtual ~IIntegration() = default;

    public:
      virtual void integrate(double from, double to, std::function<double(double)> f) = 0;

	};


  class RectangleRule : public IIntegration
  {
    public:
      RectangleRule(int N) : N(N)  {}

      void integrate(double from, double to, std::function<double(double)> f){
        double sum = 0; 
        double x, value;
        double dx = (to - from) / N;
        for(int i=1; i<=N; i++){
          x = from + i * dx;
          sum += f(x);
        }
        value = sum * dx;
        std::cout << "Integral via Rectangle Rule =  "<< value << std::endl;
      }

    private:
      int N;
  };

  class TrapezoidalRule : public IIntegration
  {
    public:
      TrapezoidalRule(int N) : N(N)  {}

      void integrate(double from, double to, std::function<double(double)> f){
        double sum = 0; 
        double a, b, value;
        double dx = (to - from) / N;
        for(int i=1; i<=N; i++){
          a = from + (i-1) * dx;
          b = from + i * dx;
          sum += (f(a)+f(b))/2;
        }
        value = sum * dx;
        std::cout << "Integral via Trapezoidal Rule =  "<< value << std::endl;
      }

    private:
      int N;
    };

  class SimpsonsRule : public IIntegration
  {
    public:
      SimpsonsRule(int N) : N(N) {}

      void integrate(double from, double to, std::function<double(double)> f){
        double value = 0, middle = 0; 
        double x;
        double dx = (to - from) / N;
        for(int i=1; i<=N; i++){
          x = from + i * dx;
          middle += f(x - dx / 2);
          if(i < N) value += f(x);
        }
        value = dx / 6 * (f(from) + f(to) + 2 * value + 4 * middle);
        std::cout << "Integral via Simpson's Rule =  "<< value << std::endl;
      }

    private:
      int N;
    };

  class MonteCarlo : public IIntegration
  {
    public:
      MonteCarlo(int samples) : samples(samples) {}

      void findPi(){
        int interval, i; 
        double x, y, d; 
        int inside = 0; 

        srand(time(NULL)); 
      
        for (int i=0; i<samples; i++) { 
      
          x = ((double) rand() / RAND_MAX);
          y = ((double) rand() / RAND_MAX);
      
          // Distance from (0,0) squared
          d = x*x + y*y; 
      
          // if true point is in the circle
          if (d <= 1) inside++;

        } 
        double pi = 4 * inside / (double)samples;
        std::cout << "[" << samples << " samples] Pi  = " << pi << std::endl;
      }

      void integrate(double from, double to, std::function<double(double)> f){

        srand(time(NULL));
        double value  = 0;
        double dx = to - from;
        for(int i=0; i<samples; i++){
          double r = from + ((double) rand()/(double)(RAND_MAX) * dx);
          value += f(r);
        }
        value = dx * value / samples;
        std::cout << "Integral via Monte Carlo method = " << value << std::endl;
      }

    private:
      int samples;
    };
}

void ex3(){

  int samples = 10;
  std::unique_ptr<integration::MonteCarlo> monteCarlo;

  for(int i=1; i<=5; i++){

      monteCarlo = std::make_unique<integration::MonteCarlo>(integration::MonteCarlo(samples));
      samples *= 10;
      monteCarlo->findPi();
    }
}


int main(void){

  auto lambda = [](double x){ return pow(x,5) + 3.5*pow(x,4) - 2.5*pow(x,3) - 12.5*pow(x,2) + 1.5*x + 9;};
  int intervals = 50;
  
  /*
  std::unique_ptr<integration::RectangleRule> rectangle = 
					std::make_unique<integration::RectangleRule>(integration::RectangleRule(intervals));

  rectangle->integrate(-3,2,lambda);

  std::unique_ptr<integration::TrapezoidalRule> trapezoidal = 
					std::make_unique<integration::TrapezoidalRule>(integration::TrapezoidalRule(intervals));

  trapezoidal->integrate(-3,2,lambda);

  std::unique_ptr<integration::SimpsonsRule> simpson = 
					std::make_unique<integration::SimpsonsRule>(integration::SimpsonsRule(intervals));

  simpson->integrate(-3,2,lambda);
  */

  std::unique_ptr<integration::MonteCarlo> monteCarlo = 
            std::make_unique<integration::MonteCarlo>(integration::MonteCarlo(intervals));

  monteCarlo->integrate(-3,2,lambda);
  
  //ex3();
  
  return 0;
}