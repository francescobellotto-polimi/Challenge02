#include "matrix.hpp"

// to speedup the implementation when T = double 

using namespace algebra;

// matrix class implementation
template class Matrix<double, StorageOrder::ROW_WISE>;
template class Matrix<double, StorageOrder::COL_WISE>;

template std::vector<double> operator * (const Matrix<double, StorageOrder::ROW_WISE> & A, 
                                               const std::vector<double> & v);
template std::vector<double> operator * (const Matrix<double, StorageOrder::COL_WISE> & A, 
                                               const std::vector<double> & v);

