#ifndef HH_MATRIX_HH
#define HH_MATRIX_HH

#include <iostream>
#include <array>
#include <map>
#include <vector>
#include <limits>
#include <optional>

namespace algebra
{

// enumerator for the storage strategy types
enum class StorageStrategy
{
    COMPRESSED,
    UNCOMPRESSED
}
;

// enumerator for the storage order types
enum class StorageOrder
{
    ROW_WISE,
    COL_WISE
}
;

// enumerator for the norm type
enum class NormType
{
    One,
    Infinity,
    Frobenius
}
;

// Map key ordering depending on the passed storage order
template <StorageOrder order=StorageOrder::ROW_WISE>
struct CompOp
{
    bool operator ()
    (const std::array<std::size_t,2> & lhs, const std::array<std::size_t,2> & rhs) const
    {
        if constexpr (order == StorageOrder::COL_WISE)
            return (lhs[1]<rhs[1] || (lhs[1]==rhs[1] && lhs[0]<rhs[0]));
        else 
            return (lhs < rhs);
    };
}
;

template <typename T, StorageOrder order=StorageOrder::ROW_WISE>
class Matrix
{
    public:
        // constructor taking number of rows and number of columns as input
        Matrix (std::size_t nrows_=0, std::size_t ncols_=0):
        m_nrows(nrows_), m_ncols(ncols_) {};

        inline bool is_compressed(void) const {return m_compressed;};

        // return the number of rows 
        inline std::size_t rows() const {return m_nrows;};

        // return the number of columns 
        inline std::size_t columns() const {return m_ncols;};

        // call operator (first parameter: row; second parameter: column)
        T & operator () (int i, int j);
        T operator () (int i, int j) const;

        // compress the matrix
        void compress(void);

        // uncompress the matrix
        void uncompress(void);

        // resize matrix
        void resize(std::size_t new_nrows, std::size_t new_ncols);

        // stream the matrix
        template <typename U, StorageOrder o> 
        friend std::ostream & operator<<(std::ostream &stream, Matrix<U, o> & M);

        // operator for the multiplication matrix - vector
        template <typename U, StorageOrder o>
        friend std::vector<U> operator * (const Matrix<U, o> & A, const std::vector<U> & v);

        // operator for the multiplication matrix - matrix (one column)
        template <typename U, StorageOrder o>
        friend std::vector<U> operator * (const Matrix<U, o> & A, const Matrix<U, o> & v);

        // compute the required norm
        template <NormType ntype = NormType::One>
        double norm() const;

        // read matrix in the matrix market format
        void read_mtx(const std::string & filename);

    private:
        // number of rows and columns
        std::size_t m_nrows;
        std::size_t m_ncols;

        // number of non zero elements
        std::size_t m_nnz=0;

        // data structure for the uncompressed storage
        std::map<std::array<std::size_t,2>, T, CompOp<order>> m_uc_data;

        // data structures for the compressed storage
        std::vector<std::size_t> m_in_idx;
        std::vector<std::size_t> m_out_idx;
        std::vector<T> m_c_values;

        // boolean indicating the storage strategy
        bool m_compressed = false;

        // method returning the index of the element we pass as parameter on m_c_values in compressed storage;
        // if not present, it returns quiet_NaN
        std::optional<std::size_t> index(std::size_t r, std::size_t c) const;

        // method updating the number of nonzero elements
        void update_nnz();
}
;

// Matrix - Vector multiplication
template <typename T, StorageOrder order>
std::vector<T> operator * (const Matrix<T, order> & A, const std::vector<T> & v);

}  // namespace algebra

// IMPLEMENTATION
#include "matrix_impl.hpp"

using namespace algebra;
// To speedup the compilation when T = double
extern template class Matrix<double, StorageOrder::ROW_WISE>;
extern template class Matrix<double, StorageOrder::COL_WISE>;
extern template std::vector<double> operator * (const Matrix<double, StorageOrder::ROW_WISE> & A, 
                                               const std::vector<double> & v);
extern template std::vector<double> operator * (const Matrix<double, StorageOrder::COL_WISE> & A, 
                                               const std::vector<double> & v);

#endif
