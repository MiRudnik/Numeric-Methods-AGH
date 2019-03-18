#include "aghMatrix.h"
#include <math.h>

// Default Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix() 
{
  rows = 0;
  cols = 0;
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(const std::vector<std::vector<T>>& mat) 
{
  matrix.resize(mat.size());
  for (unsigned i = 0; i < mat.size(); i++) 
  {
    matrix[i].resize(mat[i].size());
    for(unsigned j = 0; j < mat[i].size(); j++)
    {
      matrix[i][j] = mat[i][j];
    }
  }
  rows = matrix.size();
  cols = matrix[0].size();
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
AGHMatrix<T>::AGHMatrix(const AGHMatrix<T>& rhs) 
{
  matrix = rhs.matrix;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned AGHMatrix<T>::get_rows() const 
{
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned AGHMatrix<T>::get_cols() const 
{
  return this->cols;
}

// Assignment Operator                                                                                                                                                        
template<typename T>
AGHMatrix<T>& AGHMatrix<T>::operator=(const AGHMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  matrix.resize(new_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(new_cols);
  }

  for (unsigned i=0; i<new_rows; i++) 
  {
    for (unsigned j=0; j<new_cols; j++) 
    {
      matrix[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) 
{
  return this->matrix[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) const 
{
  return this->matrix[row][col];
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator+(const AGHMatrix<T>& rhs) 
{
  if(this->get_cols() != rhs.get_cols() || this->get_rows() != rhs.get_rows()){
    std::string txt = "Different sizes of matrices.";
    throw std::invalid_argument(txt);
  }

  AGHMatrix<T> newMatrix(this->matrix);
  for(int i=0; i<this->get_cols(); i++){
    for(int j=0; j<this->get_rows(); j++){
      newMatrix(i,j) += rhs(i,j);
    }
  }
  return newMatrix;
}

// Subtraction of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator-(const AGHMatrix<T>& rhs) 
{
  if(this->get_cols() != rhs.get_cols() || this->get_rows() != rhs.get_rows()){
    std::string txt = "Different sizes of matrices.";
    throw std::invalid_argument(txt);
  }

  AGHMatrix<T> newMatrix(this->matrix);
  for(int i=0; i<this->get_cols(); i++){
    for(int j=0; j<this->get_rows(); j++){
      newMatrix(i,j) -= rhs(i,j);
    }
  }
  return newMatrix;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator*(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement multiplication of two matrices
  if(this->get_cols() != rhs.get_rows()){
    std::string txt = "Invalid size of matrices.";
    throw std::invalid_argument(txt);
  }
  AGHMatrix<T> newMatrix(this->get_rows(),rhs.get_cols(),0);

  for(int i=0; i<this->get_rows(); i++){
    for(int j=0; j<rhs.get_cols(); j++){
      
      for(int k=0; k<this->get_cols(); k++){
        newMatrix(i,j) += this->matrix[i][k] * rhs.matrix[k][j];
      }
    }
  }
  return newMatrix;
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i=0; i<matrix.get_rows(); i++) 
  { 
    for (int j=0; j<matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
    stream << std::endl;
}

template<typename T>
bool AGHMatrix<T>::isSymmetric()
{
  if(this->get_cols() != this->get_rows()){
    return false;
  }
  for (int i=0; i<this->get_rows(); i++){ 
    for (int j=0; j<this->get_cols(); j++){
      if(this->matrix[i][j] != this->matrix[j][i]) return false;
    }
  }
  return true;
}

template<typename T>
T detRecursively(std::vector<std::vector<T>> matrix, int size)
{
  // Base cases 
  if(size == 1) return matrix[0][0];
  else if(size == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  else{
    T result = 0;
    int sign = 1;
    for(int i=0; i<size; i++){
      
      // new matrix without 0-th row and i-th column
      std::vector<std::vector<double>> tmp;
      for(int j=0; j<size; j++){
        // exluding i-th row
        if(j != 0){
          tmp.push_back(matrix[j]);
          // excluding i-th column
          tmp.back().erase(tmp.back().begin() + i);
        }
      }
      result += sign * matrix[0][i] * detRecursively(tmp, size-1);
      sign = -sign;
    }
    return result;
  }
}

template<typename T>
T AGHMatrix<T>::det()
{
  if(this->get_cols() != this->get_rows()){
    std::string txt = "Not a square matrix.";
    throw std::invalid_argument(txt);
  }

  return detRecursively(this->matrix,this->get_rows());
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::transpose()
{
  AGHMatrix<T> transposed(this->get_cols(), this->get_rows(), 0.0);
  for(int i=0; i<this->get_rows(); i++){
    for(int j=0; j<this->get_cols(); j++){
      transposed(j,i) = this->matrix[i][j];
    }
  }
  return transposed;
}


template<typename T>
std::pair<AGHMatrix<T>, AGHMatrix<T>> AGHMatrix<T>::LU()
{
  if(this->get_cols() != this->get_rows()){
    std::string txt = "Not a square matrix.";
    throw std::invalid_argument(txt);
  }
  
  std::pair<AGHMatrix<T>, AGHMatrix<T>> result(
      AGHMatrix<T>(this->get_rows(),this->get_cols(),0),
      AGHMatrix<T>(this->get_rows(),this->get_cols(),0)
  );

  for(int i=0; i<result.first.get_cols(); i++){
    result.first(i,i) = 1;
  }
  T element;

  for(int i=0; i<result.first.get_cols(); i++){
    // i-th row in second matrix
    for(int j=0; j<result.first.get_rows(); j++){
      // only above or on diagonal
      if(j >= i){
        element = this->matrix[i][j];
        for(int k=0; k<i; k++){
          element -= result.first(i,k) * result.second(k,j);
        }
        result.second(i,j) = element;
      }
    }
    // i-th column in first matrix
    for(int j=0; j<result.first.get_rows(); j++){
      // only below diagonal
      if(i < j){
        element = this->matrix[j][i];
        for(int k=0; k<i; k++){
          element -= result.first(j,k) * result.second(k,i);
        }
        result.first(j,i) = element / result.second(i,i);
      }
    }
  }
  return result;
}


template<typename T>
std::pair<AGHMatrix<T>, AGHMatrix<T>> AGHMatrix<T>::cholesky()
{
  if(this->get_cols() != this->get_rows()){
    std::string txt = "Not a square matrix.";
    throw std::invalid_argument(txt);
  }
  
  AGHMatrix<T> mat(this->get_rows(),this->get_cols(),0);

  for(int i=0; i<this->get_rows(); i++){
    // only below or on diagonal
    for(int j=0; j<=i; j++){
      T sum = 0;

      if(i == j){
        for(int k=0; k<j; k++){
          sum += pow(mat(j,k), 2);
        }
        mat(j,j) = sqrt(this->matrix[j][j] - sum); 
      }
      else{
        for(int k=0; k<j; k++){
          sum += (mat(i,k) * mat(j,k));
        }
        mat(i,j) = (this->matrix[i][j] - sum) / mat(j,j); 
      }
    }
  }
  std::pair<AGHMatrix<T>, AGHMatrix<T>> result(mat,mat.transpose());
  return result;
}


template<typename T>
AGHMatrix<T> gaussEl(AGHMatrix<T> augmented)
{
  if(augmented.get_rows() != augmented.get_cols() - 1){
    std::string txt = "Incorrect size.\nAugmented matrix has size N x N+1.";
    throw std::invalid_argument(txt);
  }

  int rows = augmented.get_rows();

  // zeroing coeffs below the leading one in the i-th column
  for(int i=0; i<rows-1; i++){
    // skipping the leading coefficient
    for(int j=i+1; j<rows; j++){
        // getting the factor zeroing in the given row
        double factor = -1*augmented(j,i)/augmented(i,i);
        for(int k=i; k<rows+1; k++){
          // operation must be row elementary so repeating
          // for each element in row
          augmented(j,k) += factor*augmented(i,k);
        }
    }
  }
  // matrix is in the correct form

  for(int i=0;i<rows;i++){
		if(augmented(i,i)==0){
      std::string txt;
      if(augmented(i,i+1)==0) txt = "Infinitely many solutions.";
      else txt = "No solution.";
      throw std::invalid_argument(txt);
    }
	}
	
  AGHMatrix<T> solution(rows, 1, 0);

  // getting the answers (starting from the last row)
  for(int i=rows-1; i>=0; i--){
    // row is in form: augmented(i,i)*solution(i,1) = augmented(i,n)
    solution(i,0) = augmented(i,rows) / augmented(i,i);
    // zeroing coeffs connected with computed variable
    // (to achieve rows: a * x = b)
    for(int j=i-1; j>=0; j--){
      // its like getting the constants on the same side of equation
      // since variable has been computed
      augmented(j,rows) -= solution(i,0) * augmented(j,i);
      // not needed in the algorithm
      augmented(j,i) = 0;
    }
  }
  
  return solution;
}

