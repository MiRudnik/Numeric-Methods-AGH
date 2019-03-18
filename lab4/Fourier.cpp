#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <complex>
#include <chrono>
#include <fstream>


typedef std::vector<std::complex<double>> FourierData;

namespace Fourier
{
	class DiscreteFourier
	{
	public:
        // construction methods
        DiscreteFourier(const FourierData data) : values(data) {}
		// methods
		FourierData dft() {
            int N = values.size();

            FourierData results;

            for (int k = 0; k < N; k++) {  // For each output element
                std::complex<double> sum = 0;
                std::complex<double> tmp;
                for (int n = 0; n < N; n++) {  // For each input element
                    double angle = 2 * M_PI * n * k / N;
                    tmp.real(values[n].real() * cos(angle) + values[n].imag() * sin(angle));
                    tmp.imag(-values[n].real() * sin(angle) + values[n].imag() * cos(angle));
                    sum = sum + tmp;
                }
                results.push_back(sum);
            }
            return results;
        }
		// private members
	private:
		FourierData values;
	};

    class FastFourier
	{
		
	public:
        // construction methods
        FastFourier(const FourierData data) : values(data) {}
		// methods
		FourierData fft() {
            if(values.size()%2 != 0){
                std::string txt = "Size must be a power of two.";
                throw std::invalid_argument(txt);
            }
            FourierData result(values.size());
            computeFFT(result, 0, 0, values.size(), 1);
            return result;
        }
        // private methods
    private:
        void computeFFT(FourierData &data, int index, int from, int size, int stride){
            if(size == 1){
                data[index] = values[from];
                // one point DFT equals given point, no need to invoke dft()
                // could be done earlier so that dft is required
            }
            else{
                computeFFT(data, index, from, size/2, 2*stride);   // even elements
                computeFFT(data, index+size/2, from+stride, size/2, 2*stride);   // odd elements
                
                for(int k=index; k<index+size/2; k++){
                    std::complex<double> even = data[k];
                    std::complex<double> odd = data[k+size/2];
                    std::complex<double> tf = exp(std::complex<double>(0,-2.*M_PI*k/size));
                    data[k] = even + tf*odd;
                    data[k+size/2] = even - tf*odd;
                }
            }
        }
		// private members
	private:
        FourierData values;
	};

}



int main()
{
	
    const FourierData test4 = 
    {
        std::complex<double>(1,0),
        std::complex<double>(2,-1),
        std::complex<double>(0,-1),
        std::complex<double>(-1,2)
    };

    srand (time(NULL));
    FourierData test512;
    for(int i=0; i<4096; i++){
        std::complex<double> r(rand()%7-3,rand()%7-3);
        test512.push_back(r);
    }

    std::unique_ptr<Fourier::DiscreteFourier> smallDFT = 
					std::make_unique<Fourier::DiscreteFourier>(Fourier::DiscreteFourier(test4));

    std::chrono::high_resolution_clock::time_point dftStart = std::chrono::high_resolution_clock::now();
    FourierData results = smallDFT->dft();
    std::chrono::high_resolution_clock::time_point dftEnd = std::chrono::high_resolution_clock::now();

    auto durationDFT = std::chrono::duration_cast<std::chrono::microseconds>(dftEnd - dftStart).count()/1000.0;

   std::cout << "DFT: " << durationDFT << "ms" << std::endl;

    // FOR TASK 1
    /*for(int i=0; i<results.size();i++){
        std::cout << results[i] << std::endl;
    }*/

    std::unique_ptr<Fourier::FastFourier> smallFFT = 
					std::make_unique<Fourier::FastFourier>(Fourier::FastFourier(test4));

    std::chrono::high_resolution_clock::time_point fftStart = std::chrono::high_resolution_clock::now();
    FourierData resultsFFT = smallFFT->fft();
    std::chrono::high_resolution_clock::time_point fftEnd = std::chrono::high_resolution_clock::now();

    auto durationFFT = std::chrono::duration_cast<std::chrono::microseconds>(fftEnd - fftStart).count()/1000.0;

    std::cout << "FFT: " << durationFFT << "ms" << std::endl;

    // FOR TASK 2
    /*for(int i=0; i<resultsFFT.size();i++){
        std::cout << resultsFFT[i] << std::endl;
    } */
    
    FourierData time_series;
    // https://datamarket.com/data/set/22p5
    std::ifstream infile("time_series_data.csv");
    std::string line;
    std::getline(infile, line); // skip first line
    for(int i=0;i<128;i++){
        std::getline(infile, line);
        std::string value = line.substr(line.find(','));
        value = value.substr(1);
        value = value.substr(0,value.find(','));
        double d = std::stod(value,NULL);
        time_series.push_back(std::complex<double>(d,0));
    }
    
    std::unique_ptr<Fourier::FastFourier> time_series_analysis = 
					std::make_unique<Fourier::FastFourier>(Fourier::FastFourier(time_series));

    FourierData analysis = time_series_analysis->fft();

    for(int i=0; i<analysis.size();i++){
        std::cout << analysis[i] << std::endl;
    }

	return 0;
}