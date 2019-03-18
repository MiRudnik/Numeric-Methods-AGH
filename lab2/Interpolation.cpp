#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <memory>
#include <algorithm>
#include <ctime>

typedef std::vector<std::array<float, 2>> PointsList2D;

/// <summary>
/// Template function used for clamping specified datasets
/// </summary>
template <typename T>
T GetIndexClamped(const std::vector<T> points, int index)
{
	if (index < 0)
		return points[0];
	else if (index >= int(points.size()))
		return points.back();
	else
		return points[index];
}

namespace exceptions
{
	/// <summary>
	/// Helper class used for throwing not implemented error
	/// </summary>
	class NotImplementedException : public std::logic_error
	{
	public:
		NotImplementedException() : std::logic_error{ "Function not yet implemented." } {}
	};
}

namespace interpolation
{
	/// <summary>
	/// Defines interface for interpolation classes
	/// </summary>
	class IInterpolation
	{
		// public construction and destruction methods
	public:
		virtual ~IInterpolation() = default;

		// public interface methods
	public:
		virtual void Interpolate1D(int pointsToInterpolate) = 0;
		virtual void Interpolate2D(int pointsToInterpolate) = 0;

	};

	/// <summary>
	/// Defines Hermite Cubic Interpolation 
	/// </summary>
	class CubicInterpolation : public IInterpolation
	{
		// public construction methods
	public:
		CubicInterpolation(const PointsList2D& points) : pointsList(points) {}

		// IInterpolation methods
	public:
		void Interpolate2D(int pointsToInterpolate) override
		{
			std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0, points_size = pointsList.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				std::array<PolynomialCoeffs, 2> coeffs;
				std::array<float, 2> A = GetIndexClamped(pointsList, index[i] - 1);
				std::array<float, 2> B = GetIndexClamped(pointsList, index[i]);
				std::array<float, 2> C = GetIndexClamped(pointsList, index[i] + 1);
				std::array<float, 2> D = GetIndexClamped(pointsList, index[i] + 2);

				for (int j = 0; j < 2; j++)
				{
					coeffs[j].A = A[j];
					coeffs[j].B = B[j];
					coeffs[j].C = C[j];
					coeffs[j].D = D[j];
				}

				float x = CubicHermite(coeffs[0], t[i]);
				float y = CubicHermite(coeffs[1], t[i]);

				std::cout << "Value at " << tx[i] << " = " << x << "  " << y << std::endl;
			}
		}

		void Interpolate1D(int pointsToInterpolate) override
		{
			// TASK1
			std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0, points_size = pointsList.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				PolynomialCoeffs coeffs;
				std::array<float, 2> A = GetIndexClamped(pointsList, index[i] - 1);
				std::array<float, 2> B = GetIndexClamped(pointsList, index[i]);
				std::array<float, 2> C = GetIndexClamped(pointsList, index[i] + 1);
				std::array<float, 2> D = GetIndexClamped(pointsList, index[i] + 2);

				
				coeffs.A = A[0];
				coeffs.B = B[0];
				coeffs.C = C[0];
				coeffs.D = D[0];
				

				float x = CubicHermite(coeffs, t[i]);

				std::cout << "Value at " << tx[i] << " = " << x << std::endl;
			}
		}

		// private methods
	private:
		struct PolynomialCoeffs
		{
			float A, B, C, D, t;
		};

		float CubicHermite(PolynomialCoeffs coeffs, float t) const
		{
			float a = -coeffs.A / 2.0f + (3.0f*coeffs.B) / 2.0f - (3.0f*coeffs.C) / 2.0f + coeffs.D / 2.0f;
			float b = coeffs.A - (5.0f*coeffs.B) / 2.0f + 2.0f*coeffs.C - coeffs.D / 2.0f;
			float c = -coeffs.A / 2.0f + coeffs.C / 2.0f;
			float d = coeffs.B;

			return a * pow(t, 3) + b * pow(t, 2) + c * t + d;
		}

		// private members
	private:
		const PointsList2D& pointsList;

	};

	class LagrangeInterpolation : public IInterpolation
	{
	public:
		LagrangeInterpolation(const PointsList2D& points) : pointsList(points) {}

		void Interpolate1D(int pointsToInterpolate) override
		{
			// TASK2
			std::vector<int> index(pointsToInterpolate);
			std::vector<float> x;

			int i = 0, points_size = pointsList.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &x]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				x.push_back((points_size)* percent);
				return int(x[i++]);
			});

			for(int k = 0; k < pointsToInterpolate; k++){
				// loop for each value to calculate
				float result = 0.0f;
				
				for (int i = 0; i < points_size; i++)
				{
					// summing

					float product = pointsList[i][2];
					for (int j = 0; j < points_size; j++)
					{
						// multiplicating
						if (j!=i) product *= (x[k] - pointsList[j][0]) / double(pointsList[i][0] - pointsList[j][0]);
					}

					result += product;
					
				}

				std::cout << "Value at " << x[k] << " = " << result << std::endl;
			}
			
		};
		void Interpolate2D(int pointsToInterpolate) override
		{
			// TASK2 : TODO
			// LUL NOPE
		};

		// private members
	private:
		const PointsList2D& pointsList;
	};
}

namespace regression
{
	// Mean Squared Error
	float mse(std::vector<float> x, std::vector<float> y)
	{
		// TASK3
		if(x.size() != y.size()){
			return NAN;
		}
		
		float result = 0.0f;
		for(int i=0;i<x.size();i++){
			result += pow((x[i]-y[i]),2);
		}
		return result/x.size();
	}

	std::pair<float, float> linearRegression(std::vector<float> x, std::vector<float> y)
	{
		// TASK4
		std::pair<float, float> coeffs;	// y = ax + b
		if(x.size() != y.size()){
			coeffs.first = NAN;
			coeffs.second = NAN;
			return coeffs;
		}

		int size = x.size();
		float meanx=0, meany=0;
		for(int i=0;i<size;i++){
			meanx += x[i];
			meany += y[i];
		}
		meanx /= size;
		meany /= size;

		std::vector<float> errorx (size);
		std::vector<float> errory (size);
		float num = 0.0f;	// sum(errorx*errory)
		float denom = 0.0f;	// sum(errorx^2)
		for(int i=0;i<size;i++){
			errorx[i] = x[i] - meanx;
			errory[i] = y[i] - meany;
			num += errorx[i]*errory[i];
			denom += pow(errorx[i],2);
		}

		// a = sum(errorx * errory) / sum(errorx^2)
		coeffs.first = num/denom;
		coeffs.second = meany - coeffs.first * meanx;

		return coeffs;
	}

	class LinearRegressor
	{
	public:
		// TASK6 : TODO
		LinearRegressor() { }
		float fit(std::vector<float> x, std::vector<float> y){ throw exceptions::NotImplementedException(); }
		std::vector<float> predict(std::vector<float> x) { throw exceptions::NotImplementedException(); }

	private:
		std::pair<float, float> coefficients;
	};
}

int main()
{
	
	const PointsList2D points2D =
	{
		{ 0.0f, 1.1f },
		{ 1.6f, 8.3f },
		{ 2.3f, 6.5f },
		{ 3.5f, 4.7f },
		{ 4.3f, 3.1f },
		{ 5.9f, 7.5f },
		{ 6.8f, 0.0f }
	};

	std::unique_ptr<interpolation::CubicInterpolation> interpolationH = 
					std::make_unique<interpolation::CubicInterpolation>(interpolation::CubicInterpolation(points2D));

	std::cout << std::endl << "Hermit Interpolation 2D" << std::endl;
	interpolationH->Interpolate2D(13);

	std::cout << std::endl << "Hermit Interpolation 1D" << std::endl;
	interpolationH->Interpolate1D(13);

	std::unique_ptr<interpolation::LagrangeInterpolation> interpolationL = 
					std::make_unique<interpolation::LagrangeInterpolation>(interpolation::LagrangeInterpolation(points2D));

	std::cout << std::endl << "Lagrange Interpolation 1D" << std::endl;
	interpolationL->Interpolate1D(13);

	std::cout << std::endl << "Mean Squared Error" << std::endl;
	std::vector<float> v1 (10);
	std::vector<float> v2 (10);
	srand((time(0)));
	std::generate (v1.begin(), v1.end(), [](){return std::rand()/(RAND_MAX/10.0);});
	std::generate (v2.begin(), v2.end(), [](){return std::rand()/(RAND_MAX/10.0);});
	for(int i=0;i<v1.size();i++){
		std::cout << v1[i] << "\t" << v2[i] << std::endl;
	}
	std::cout << std::endl << "MSE = " << regression::mse(v1,v2) << std::endl << std::endl;

	std::vector<float> v3 {10,20,40,45,60,65,75,80};
	std::vector<float> v4 {40,40,60,80,90,110,100,130};
	std::pair<float, float> coeffs (regression::linearRegression(v3,v4));
	std::cout << "Linear regression" << std::endl;
	std::cout <<  "y = " << coeffs.first << "*x + " << coeffs.second << std::endl;


	return 0;
}