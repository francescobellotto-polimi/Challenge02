#ifndef HH_MATRIX_IMPL_HH
#define HH_MATRIX_IMPL_HH

#include <fstream>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <utility>
#include <cmath>
#include "matrix.hpp"

namespace algebra

{

template <typename T, StorageOrder order>
void Matrix<T, order>::update_nnz()
{
    if (!m_compressed) m_nnz = m_uc_data.size();
    else m_nnz = m_c_values.size();
}


template <typename T, StorageOrder order>
std::optional<std::size_t> Matrix<T, order>::index(std::size_t r, std::size_t c) const
{
    if constexpr (order == StorageOrder::ROW_WISE)
    {
        if (m_in_idx[r]!=m_in_idx[r+1] && 
        (std::find(m_out_idx.begin()+m_in_idx[r], m_out_idx.begin()+m_in_idx[r+1],c)!=m_out_idx.begin()+m_in_idx[r+1]))

            return m_in_idx[r] + 
            std::distance(m_out_idx.begin(),std::find(m_out_idx.begin()+m_in_idx[r], m_out_idx.begin()+m_in_idx[r+1],c));
    }

    else
    {
        if (m_in_idx[c]!=m_in_idx[c+1] && 
        (std::find(m_out_idx.begin()+m_in_idx[c], m_out_idx.begin()+m_in_idx[c+1],c)!=m_out_idx.begin()+m_in_idx[c+1]))

            return m_in_idx[c] + 
            std::distance(m_out_idx.begin(),std::find(m_out_idx.begin()+m_in_idx[c], m_out_idx.begin()+m_in_idx[c+1],r));
    }

    // create an uninitialized optional size_t, so that the call to index produces a "boolean" set to false
    std::optional<std::size_t> uninit;
    return uninit;
}
;


template <typename T, StorageOrder order>
T & Matrix<T, order>::operator () (int i, int j)
{
    std::size_t r = static_cast<std::size_t>(i);
    std::size_t c = static_cast<std::size_t>(j);

    if (!m_compressed)
    {
        auto it = m_uc_data.find({r,c});
        if (it==m_uc_data.cend())
        {
            m_uc_data[{r,c}] = 0;
            return m_uc_data[{r,c}];
        }

        return it->second;
    }

    else
    {
        auto idx = index(r,c);
        if (idx)
            return m_c_values[*idx];

        std::cerr << "Non-const call operator was asked for uninitialized element: decompress matrix first" << std::endl;
        std::exit(1);
    }
}
;


template <typename T, StorageOrder order>
T Matrix<T, order>::operator () (int i, int j) const
{
    std::size_t r = static_cast<std::size_t>(i);
    std::size_t c = static_cast<std::size_t>(j);

    if (!m_compressed)
    {
        const auto it = m_uc_data.find({r,c});
        if (it==m_uc_data.cend())
            return static_cast<T>(0);
        return it->second;
    }

    else
    {
        auto idx = index(r,c);
        if (idx)
            return m_c_values[*idx];

        return static_cast<T>(0);
    }
}
;

template <typename T, StorageOrder order>
void Matrix<T, order>::compress(void)
{
    update_nnz();

    // if the matrix is already compressed the function does not do anything
    if (m_compressed) return;

    // reserve space for the new data structures
    m_in_idx.reserve(m_nrows+1);
    m_out_idx.reserve(m_uc_data.size());
    m_c_values.reserve(m_uc_data.size());

    if constexpr (order == StorageOrder::ROW_WISE)
    {
        std::size_t k = 0;

        for (std::size_t r=0; r<m_nrows; ++r)
        {
            // find the index from which row r starts
            m_in_idx.emplace_back(k);

            for (auto it = m_uc_data.lower_bound({r,0}); it != m_uc_data.upper_bound({r,m_ncols-1}); ++it)
            {
                // fill the out index vector with column indexes related to row r
                m_out_idx.emplace_back(it->first[1]);

                // fill with the values related to row r
                m_c_values.emplace_back(it->second);

                ++k;
            }
        }

        m_in_idx.emplace_back(k);
    }

    else
    {
        std::size_t k = 0;
        
        for (std::size_t c=0; c<m_ncols; ++c)
        {
            // find the index from which column c starts
            m_in_idx.emplace_back(k);

            for (auto it = m_uc_data.lower_bound({0,c}); it != m_uc_data.upper_bound({m_nrows-1,c}); ++it)
            {
                // fill the out index vector with row indexes related to column c
                m_out_idx.emplace_back(it->first[0]);

                // fill with the values related to column c
                m_c_values.emplace_back(it->second);

                ++k;
            }
        }

        m_in_idx.emplace_back(k);
    }

    // clear the old data structures
    m_uc_data.clear();

    // change state of the object
    m_compressed = true;

}

template <typename T, StorageOrder order>
void Matrix<T, order>::uncompress(void)
{
    update_nnz();

    // if the matrix is already uncompressed the function does not do anything
    if (!m_compressed) return;

    if constexpr (order == StorageOrder::ROW_WISE)
    {
        for (std::size_t r=0; r<m_nrows; ++r)
        {
            if (m_in_idx[r] != m_in_idx[r+1])
            {
                for (std::size_t j=m_in_idx[r]; j<m_in_idx[r+1]; ++j)
                {
                    std::size_t c = m_out_idx[j];
                    m_uc_data[{r,c}] = m_c_values[j];
                }
            }
        }
    }

    else
    {
        for (std::size_t c=0; c<m_ncols; ++c)
        {
            if (m_in_idx[c] != m_in_idx[c+1])
            {
                for (std::size_t i=m_in_idx[c]; i<m_in_idx[c+1]; ++i)
                {
                    std::size_t r = m_out_idx[i];
                    m_uc_data[{r,c}] = m_c_values[i];
                }
            }
        }
    }

    // free the old data structures
    m_in_idx.clear();
    m_out_idx.clear();
    m_c_values.clear();

    // change state of the object
    m_compressed = false;

}

template <typename T, StorageOrder order>
void Matrix<T, order>::resize(std::size_t new_nrows, std::size_t new_ncols)
{
    uncompress();

    if (m_nrows<=new_nrows && m_ncols<=new_ncols)
    {
        m_nrows = new_nrows;
        m_ncols = new_ncols;
    }

    else
    {
        std::erase_if (m_uc_data, [new_nrows, new_ncols](const auto & item)
        {
            auto const & [key, value] = item;
            return (key[0]>=new_nrows || key[1]>=new_ncols);
        });
        m_nrows = new_nrows;
        m_ncols = new_ncols;
    }

    update_nnz();
}

template <typename T, StorageOrder order>
std::ostream & operator << (std::ostream &stream, Matrix<T, order> & M)
{
    M.update_nnz();
    stream << "Rows: " << M.m_nrows << "; " << 
              "Columns: " << M.m_ncols << "; " <<
              "Number of non-zeros: " << M.m_nnz << std::endl;

    stream << "Matrix: [" << std::endl;

    if (!M.m_compressed)
    {
        for (auto it = M.m_uc_data.cbegin(); it != M.m_uc_data.cend(); ++it)
        {
            stream << "(" << it->first[0]+1 << 
                      "," << it->first[1]+1 << 
                      "): " << it->second << std::endl;
        }
    }

    else
    {
        if constexpr (order == StorageOrder::ROW_WISE)
        {
            for (std::size_t r=0; r<M.m_nrows; ++r)
            {
                if (M.m_in_idx[r]!=M.m_in_idx[r+1])
                {
                    for (auto pos=M.m_in_idx[r]; pos<M.m_in_idx[r+1]; ++pos)
                    {
                        stream << "(" << r+1 << 
                                  "," << M.m_out_idx[pos]+1 << 
                                  "): " << M.m_c_values[pos] << std::endl;
                    }
                }
            }
        }

        else
        {
            for (std::size_t c=0; c<M.m_ncols; ++c)
            {
                if (M.m_in_idx[c]!=M.m_in_idx[c+1])
                {
                    for (auto pos=M.m_in_idx[c]; pos<M.m_in_idx[c+1]; ++pos)
                    {
                        stream << "(" << M.m_out_idx[pos]+1 << 
                                  "," << c+1 << 
                                  "): " << M.m_c_values[pos] << std::endl;
                    }
                }
            }
        }
    }

    stream << "]" << std::endl;

    return stream;
}

template <typename T, StorageOrder order>
std::vector<T> operator * (const Matrix<T, order> & A, const std::vector<T> & v)
{
    // abort if matrix size are inconsistent
    if (A.m_ncols != v.size())
    {
        std::cerr << "Inconsistent size for matrix-vector multiplication" << std::endl;
        std::exit(2);
    }

    // creation of the result
    std::vector<T> b(A.m_nrows);
    for (std::size_t r=0; r<A.m_nrows; ++r)
    {
        b[r] = static_cast<T>(0);
    }

    // case row-wise storage
    if constexpr (order==StorageOrder::ROW_WISE)
    {
        if (!A.m_compressed)
        {
            for (std::size_t r=0; r<A.m_nrows; ++r)
            {
                for (auto it = A.m_uc_data.lower_bound({r,0}); it !=  A.m_uc_data.upper_bound({r,A.m_ncols-1}); ++it)
                    b[r] += it->second * v[it->first[1]];
            }
        }

        else
        {
            for (std::size_t r=0; r<A.m_nrows; ++r)
            {
                for (std::size_t j=A.m_in_idx[r]; j<A.m_in_idx[r+1]; ++j)
                {
                    b[r] += A.m_c_values[j] * v[A.m_out_idx[j]];
                }
            }
        }
    }

    // case column-wise storage
    else 
    {
        if (!A.m_compressed)
        {
            for (std::size_t c=0; c<A.m_ncols; ++c)
            {
                for (auto it = A.m_uc_data.lower_bound({0,c}); it !=  A.m_uc_data.upper_bound({A.m_nrows-1,c}); ++it)
                {
                    b[it->first[0]] += it->second * v[c]; 
                }
            }
        }

        else
        {
            for (std::size_t c=0; c<A.m_ncols; ++c)
            {
                for (std::size_t i=A.m_in_idx[c]; i<A.m_in_idx[c+1]; ++i)
                {
                    b[A.m_out_idx[i]] += A.m_c_values[i] * v[c];
                }
            }
        }
    }

    return b;
}


template <typename T, StorageOrder order>
std::vector<T> operator * (const Matrix<T, order> & A, const Matrix<T, order> & v)
{
    // abort if matrix size are inconsistent or the second matrix is not a column
    if (A.m_ncols != v.m_nrows || v.m_ncols!=1)
    {
        std::cerr << "Inconsistent size for matrix-vector multiplication" << std::endl;
        std::exit(2);
    }

    // creation of the result
    std::vector<T> b(A.m_nrows);
    for (std::size_t r=0; r<A.m_nrows; ++r)
    {
        b[r] = static_cast<T>(0);
    }
    
    // case row-wise storage
    if constexpr (order==StorageOrder::ROW_WISE)
    {
        if (!A.m_compressed)
        {
            for (std::size_t r=0; r<A.m_nrows; ++r)
            {
                for (auto it = A.m_uc_data.lower_bound({r,0}); it !=  A.m_uc_data.upper_bound({r,A.m_ncols-1}); ++it)
                    b[r] += it->second * v(it->first[1],0);
            }
        }

        else
        {
            for (std::size_t r=0; r<A.m_nrows; ++r)
            {
                for (std::size_t j=A.m_in_idx[r]; j<A.m_in_idx[r+1]; ++j)
                {
                    b[r] += A.m_c_values[j] * v(A.m_out_idx[j],0);
                }
            }
        }
    }

    // case column-wise storage
    else 
    {
        if (!A.m_compressed)
        {
            for (std::size_t c=0; c<A.m_ncols; ++c)
            {
                for (auto it = A.m_uc_data.lower_bound({0,c}); it !=  A.m_uc_data.upper_bound({A.m_nrows-1,c}); ++it)
                {
                    b[it->first[0]] += it->second * v(c,0); 
                }
            }
        }

        else
        {
            for (std::size_t c=0; c<A.m_ncols; ++c)
            {
                for (std::size_t i=A.m_in_idx[c]; i<A.m_in_idx[c+1]; ++i)
                {
                    b[A.m_out_idx[i]] += A.m_c_values[i] * v(c,0);
                }
            }
        }
    }

    return b;    
}


template <typename T, StorageOrder order>
template <NormType ntype>
double Matrix<T, order>::norm() const
{
    if constexpr (ntype == NormType::One)
    {
        // define the operation to sum the absolute value of elements in a map
        auto sum_abs = [] (double value, const auto & p){return std::abs(value) + static_cast<double>(std::abs(p.second));};

        // vector where we will save the sum of the elements in each column
        std::vector<double> col_sums(m_ncols,static_cast<double>(0));

        if constexpr (order == StorageOrder::ROW_WISE)
        {
            if (!m_compressed)
            {
                // loop over all elements of the map and sum the value in the right position of col_sums
                for (const auto & val: m_uc_data)
                {
                    col_sums[val.first[1]] += std::abs(val.second);
                }
            }

            else
            {
                for (std::size_t i=0; i<m_out_idx.size(); ++i)
                {
                    col_sums[m_out_idx[i]] += std::abs(m_c_values[i]);
                }   
            }
        }

        else
        {
            if (!m_compressed)
            {
                // loop over all columns to sum 
                for (std::size_t c=0; c<m_ncols; ++c)
                {   
                    // compute the sum for column c
                    col_sums[c] = std::accumulate(m_uc_data.lower_bound({0,c}), 
                                                  m_uc_data.upper_bound({m_nrows-1,c}),
                                                  0.,
                                                  sum_abs);
                }
            }

            else
            {
                for (std::size_t c=0; c<m_ncols; ++c)
                {
                    for (std::size_t i=m_in_idx[c]; i<m_in_idx[c+1]; ++i)
                    {
                        col_sums[c] += std::abs(m_c_values[i]);
                    }
                }
            }
        }

        // get the maximum among all the columns
        return std::ranges::max(col_sums);
    }


    if constexpr (ntype == NormType::Infinity)
    {
        // define the operation to sum the absolute value of elements in a map
        auto sum_abs = [] (double value, const auto & p){return std::abs(value) + static_cast<double>(std::abs(p.second));};

        // vector where we will save the sum of the elements in each row
        std::vector<double> row_sums(m_nrows,static_cast<double>(0));

        if constexpr (order == StorageOrder::ROW_WISE)
        {
            if (!m_compressed)
            {
                // loop over all rows to sum 
                for (std::size_t r=0; r<m_nrows; ++r)
                {   
                    // compute the sum for row r
                    row_sums[r] = std::accumulate(m_uc_data.lower_bound({r,0}), 
                                                  m_uc_data.upper_bound({r,m_ncols-1}),
                                                  0.,
                                                  sum_abs);
                }
            }

            else
            {
                for (std::size_t r=0; r<m_nrows; ++r)
                {
                    for (std::size_t i=m_in_idx[r]; i<m_in_idx[r+1]; ++i)
                    {
                        row_sums[r] += std::abs(m_c_values[i]);
                    }
                }
            }
        }

        else
        {
            if (!m_compressed)
            {
                // loop over all elements of the map and sum the value in the right position of row_sums
                for (const auto & val: m_uc_data)
                {
                    row_sums[val.first[0]] += std::abs(val.second);
                }
            }

            else
            {
                for (std::size_t i=0; i<m_out_idx.size(); ++i)
                {
                    row_sums[m_out_idx[i]] += std::abs(m_c_values[i]);
                }   
            }
        }

        // get the maximum among all the rows
        return std::ranges::max(row_sums);
    }


    if constexpr (ntype == NormType::Frobenius)
    {
        if (!m_compressed)
        {
            // define the operation to sum the absolute value of elements in a map
            auto sum_square = [] (double value, const auto & p){return std::abs(value) + (static_cast<double>(std::abs(p.second))
                                                                                   * static_cast<double>(std::abs(p.second)));};

            return std::sqrt(std::accumulate(m_uc_data.cbegin(), 
                             m_uc_data.cend(),
                             0.,
                             sum_square));
        }

        else
        {
            // define the operation to sum the absolute value of elements in a vector
            auto sum_square = [] (double value, const auto & p){return std::abs(value) + (static_cast<double>(std::abs(p))
                                                                                   * static_cast<double>(std::abs(p)));};

            return std::sqrt(std::accumulate(m_c_values.cbegin(), 
                             m_c_values.cend(),
                             0.,
                             sum_square));
        }
    }
}


template <typename T, StorageOrder order>
void Matrix<T, order>::read_mtx(const std::string & filename)
{
    std::ifstream file(filename);

    std::size_t nrows, ncols, nlines;

    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> nrows >> ncols >> nlines;

    // Uncompress matrix, clear it and resize it
    uncompress();
    resize(0,0);
    resize(nrows,ncols);

    T data;
    std::size_t row,col;

    // Fill the matrix 
    for (std::size_t l=0; l<nlines; ++l)
    {
        file >> row >> col >> data;
        (*this)(row-1, col-1) = data;
    }

    update_nnz();

    file.close();
}

}

#endif
