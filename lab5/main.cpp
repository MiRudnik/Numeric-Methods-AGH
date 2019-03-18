#include "aghMatrix.h"
#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>

namespace LinearEquations
{
  class EquationSolver
  {
    public:
      EquationSolver(AGHMatrix<double> matrix, AGHMatrix<double> rightSide) : b(rightSide){
        if(matrix.get_rows() != matrix.get_cols() || matrix.get_cols() != rightSide.get_rows()){
          std::string txt = "Invalid sizes.";
          throw std::invalid_argument(txt);
        }
        decompose(matrix);
      }

      AGHMatrix<double> Jacobi(int iter, AGHMatrix<double> x){
        
        AGHMatrix<double> newX = x;
        // for each x 
        for(int i=0; i<x.get_rows(); i++){
          double sum = 0.0;
          // summing
          for(int j=0; j<L.get_cols(); j++){
            if(i>j) sum += L(i,j)*x(j,0);
            else if(i<j) sum += U(i,j)*x(j,0);
          }
          newX(i,0) = (b(i,0) - sum)/D(i,i);
        }
        
        std::cout << std::setprecision(10) << "Iteration: " << iter << std::endl;
        for(int i=0; i<x.get_rows(); i++){
          std::cout << newX(i,0) << ", ";
        }
        std::cout << std::endl;

        bool check = true;
        for(int i=0; i<newX.get_rows(); i++){
          if(abs(newX(i,0) - x(i,0)) > 1e-10){
            check = false;
            break;
          }
        }
        if(check == true) return newX;
        else Jacobi(iter+1,newX);
      }

      AGHMatrix<double> GaussSeidel(int iter, AGHMatrix<double> x){
        
        AGHMatrix<double> newX = x;
        // for each x 
        for(int i=0; i<x.get_rows(); i++){
          double sum = 0.0;
          // summing
          for(int j=0; j<L.get_cols(); j++){
            // newX in i-th iteration has new elements in 0:i-1 and old ones in i:N-1
            if(i>j) sum += L(i,j)*newX(j,0);
            else if(i<j) sum += U(i,j)*newX(j,0);
          }
          newX(i,0) = (b(i,0) - sum)/D(i,i);
        }
        
        std::cout << std::setprecision(10) << "Iteration: " << iter << std::endl;
        for(int i=0; i<x.get_rows(); i++){
          std::cout << newX(i,0) << ", ";
        }
        std::cout << std::endl;

        bool check = true;
        for(int i=0; i<newX.get_rows(); i++){
          if(abs(newX(i,0) - x(i,0)) > 1e-10){
            check = false;
            break;
          }
        }
        if(check == true) return newX;
        else GaussSeidel(iter+1,newX);
      }

      AGHMatrix<double> SOR(int iter, AGHMatrix<double> x, double relaxation){
        
        AGHMatrix<double> newX = x;
        // for each x 
        for(int i=0; i<x.get_rows(); i++){
          double sum = 0.0;
          // summing
          for(int j=0; j<L.get_cols(); j++){
            // newX in i-th iteration has new elements in 0:i-1 and old ones in i:N-1
            if(i>j) sum += L(i,j)*newX(j,0);
            else if(i<j) sum += U(i,j)*newX(j,0);
          }
          newX(i,0) = (1-relaxation)*newX(i,0) + (b(i,0) - sum)/D(i,i)*relaxation;
        }
        
        std::cout << "Iteration: " << iter << std::endl;
        for(int i=0; i<x.get_rows(); i++){
          std::cout << std::setprecision(10) << newX(i,0) << ", ";
        }
        std::cout << std::endl;

        bool check = true;
        for(int i=0; i<newX.get_rows(); i++){
          if(abs(newX(i,0) - x(i,0)) > 1e-10){
            check = false;
            break;
          }
        }
        if(check == true) return newX;
        else SOR(iter+1,newX, relaxation);
      }

    private:
      AGHMatrix<double> L;
      AGHMatrix<double> D;
      AGHMatrix<double> U;
      AGHMatrix<double> b;

      void decompose(AGHMatrix<double> A){
        
        AGHMatrix<double> lower(A.get_cols(), A.get_cols(), 0.0);
        AGHMatrix<double> diagonal(A.get_cols(), A.get_cols(), 0.0);
        AGHMatrix<double> upper(A.get_cols(), A.get_cols(), 0.0);

        for(int i=0; i<A.get_rows(); i++){
          for(int j=0; j<A.get_cols(); j++){
            if(i>j) lower(i,j) = A(i,j);
            else if(i<j) upper(i,j) = A(i,j);
            else diagonal(i,j) = A(i,j);
          }
        }
        L = lower;
        D = diagonal;
        U = upper;
      }
  };
} // namespace LinearEquations

int main()
{
  std::vector<std::vector<double>> init{{4.0, 1.0, 3.0},
                                        {1.0, 5.0, 1.0},
                                        {2.0, -1.0, 8.0}};

  AGHMatrix<double> mat1A(init);

  AGHMatrix<double> mat1b({{17.0, 14.0, 12.0}});
  mat1b = mat1b.transpose();

  AGHMatrix<double> mat1x({{0.0, 0.0, 0.0}});
  mat1x = mat1x.transpose();

  std::unique_ptr<LinearEquations::EquationSolver> m1Solver = 
					std::make_unique<LinearEquations::EquationSolver>(LinearEquations::EquationSolver(mat1A, mat1b));

  std::cout << "Jacobi" << std::endl;
  AGHMatrix<double> m1res1 = m1Solver->Jacobi(1,mat1x);
  std::cout << "Gauss-Seidel" << std::endl;
  AGHMatrix<double> m1res2 = m1Solver->GaussSeidel(1,mat1x);
  std::cout << "SOR" << std::endl;
  AGHMatrix<double> m1res3 = m1Solver->SOR(1,mat1x,1.06);
  
  return 0;
}