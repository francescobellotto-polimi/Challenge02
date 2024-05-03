# IMPLEMENTATION OF SPARSE MATRIX #

The code implements a sparse matrix data structure throughout a template class which is called Matrix.

- The two template parameters are: type of the data stored in the matrix (it can be either an integral type or std::complex), the storage order, which can be either ROW_WISE or COL_WISE. The difference among the two storage strategies can be only seen when the matrix is output through the stream operator (<<): in the former case, the matrix will be traversed row-wise, in the latter column-wise.

- Matrix class DOES NOT implement range checking. The number of rows and columns may be specified when the constructor is called or through the method "resize", but one could still access elemnts out of range without getting any warning (but of course this might cause unexpected behaviour).

- Objects of class Matrix can be in one of two states: compressed and uncompressed. When created through the constructor, the matrix is uncompressed, but it can be compressed and uncompressed thanks to the two methods "compress" and "uncompress". In the compressed state, matrix-vector multiplication is faster, but the non-const call operator cannot be used to modify the matrix (i.e., if A is an object of the Matrix class in compressed state, one can not call A(i,j) if position (i,j) has not been already filled while A was in the uncompressed state): only the const version can be called.

- Another way of initializing objects of class Matrix is to read from input file a matrix in the matrix market format; this can be done through the method "read_mtx", which accepts as input the name of a file with ".mtx" extension.

- One can multiply row-by-column objects of the Matrix class simply using the "*" operator, putting at the left-hand-side a matrix and at the right-hand-side an std::vector of data of the same type T, or another object of class Matrix, provided that it only has one column, the data it stores is of type T and it is stored according to the same order of the matrix. Of course, the two dimensions must be compatible with the row-by-column product

- Finally, a template method "norm" is implemented, which accepts as template parameter the type of norm we want to compute: One, Infinity or Frobenius (defined in the "algebra" namespace, as enumerator class "NormType").

- The file matrix_double.cpp is useful just to externally implement the Matrix class in the particular case T = double (to speedup the code in such a case).


# HOW TO RUN THE CODE #

The code can be simply run by doing "make": this will create an executable which is called "main". 
In the Makefile, variable "PACS_ROOT" has to be changed with the path where the folder pacs-examples is situated (the reason why we need to access it is to use "chrono" utility to measure the speed of the operations).
In the same folder of the main.cpp, there must be a file called "lnsp_131.mtx", in order to test the implementation on the lnsp_131 matrix (to be stored in matrix market format).

# THE MAIN #

In the main.cpp file, we first create a matrix with ROW_WISE storage and one with COL_WISE storage, and initialize them through the "Insp_131.mtx" matrix (that can be found at https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html). We create a vector of doubles that are randomly generated, and then perform matrix-vector multiplication for different storage strategies (compressed and uncompressed). Then we display the time it took to perform the multiplication for all the four cases.

In the second part, we compute the three different norms according to the different storage strategies, and once again print out the time it took to perform all the required computations.

Finally, there is a commented part where everything is done once again but storing complex data (in the matrix we simply store the same real values as before but in a std::complex fashion, while in the vector we really put random complex numbers).
