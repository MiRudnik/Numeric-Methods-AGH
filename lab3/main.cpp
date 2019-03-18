#include "aghMatrix.h"
#include <iostream>

int main() 
{
    // initialize matrices using init value
    AGHMatrix<double> mat1(5, 5, 1.2);
    AGHMatrix<double> mat2(5, 5, 2.8);

    AGHMatrix<double> mat1plus2 = mat1 + mat2;
    std::cout << "Addition:" << std::endl << mat1plus2;

    AGHMatrix<double> mat3({{1,2,3},
                          {4,5,6}});
    AGHMatrix<double> mat4({{7,8},
                          {9,10},
                          {11,12}});
    
    AGHMatrix<double> mat3mul4 = mat3 * mat4;
    std::cout << "Multiplication:" << std::endl << mat3mul4;

    // initialize matrix using specified values
    std::vector<std::vector<double>> init { { 1.0, 2.0, 3.0 }, 
                                            { 2.0, 5.0, 4.0 }, 
                                            { 3.0, 4.0, 9.0 } }; 

    AGHMatrix<double> mat5(init);
    std::cout << "Matrix:" << std::endl << mat5;
    std::cout << "Is symmetric? " << mat5.isSymmetric() << std::endl;
    std::cout << "Det = " << mat5.det() << std::endl;

    AGHMatrix<double> tmat3 = mat3.transpose();
    std::cout << std::endl << "Transposition: " << std::endl;
    std::cout << "Before:" << std::endl << mat3;
    std::cout << "After:" << std::endl << tmat3;


    std::vector<std::vector<double>> init_LU{{ 5.0, 3.0, 2.0 }, 
                                            { 1.0, 2.0, 0.0 }, 
                                            { 3.0, 0.0, 4.0 }};
    AGHMatrix<double> mat6(init_LU);

    std::pair<AGHMatrix<double>, AGHMatrix<double>> lu = mat6.LU();
    std::cout << "LU decomposition: " << std::endl << lu.first << lu.second;


    std::vector<std::vector<double>> init_cholesky {{ 4.0, 12.0, -16.0 }, 
                                                    { 12.0, 37.0, -43.0 }, 
                                                    { -16.0, -43.0, 98.0 }};
    AGHMatrix<double> mat7(init_cholesky);

    std::pair<AGHMatrix<double>, AGHMatrix<double>> cholesky = mat7.cholesky();
    std::cout << "Cholesky decomposition: " << std::endl << cholesky.first << cholesky.second;

    std::vector<std::vector<double>> init_gauss {{0.0001, -5.0300, 5.8090, 7.8320, 9.5740},
                                                 {2.2660, 1.9950,  1.2120, 8.0080, 7.2190},
                                                 {8.8500, 5.6810,  4.5520, 1.3020, 5.7300},
                                                 {6.7750, -2.253,  2.9080, 3.9700, 6.2910}};
    AGHMatrix<double> augmented(init_gauss);

    AGHMatrix<double> solution = gaussEl(augmented);

    std::cout << "Gauss solution: "<< std::endl << solution;

    /*  Result
        [[ 0.21602477]
        [-0.00791511]
        [ 0.63524333]
        [ 0.74617428]]*/


    return 0;
}