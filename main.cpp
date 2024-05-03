#include <iostream>
#include <complex>
#include <random>

#include "matrix.hpp"
#include "chrono.hpp"


int main(int argc, char**argv)
{
    Timings::Chrono clock;
    std::string filename = "lnsp_131.mtx";



    /////////////////////////////// DATA ARE DOUBLES ///////////////////////////////

    std::cout << "DOUBLE DATA\n\n\n"  << std::endl; 
    
    using data_type = double;

    // reads the required matrix
    algebra::Matrix<data_type> A;
    A.read_mtx(filename);

    // generates the vector to multiply
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    std::vector<data_type> v(A.columns());
    for (std::size_t i=0; i<A.columns(); ++i)
        v[i] = uniform(generator);
    



    /////////////////////////////// TEST THE PRODUCT MATRIX - VECTOR ///////////////////////////////
    std::cout << "\n\nTESTING MATRIX - VECTOR PRODUCT" << std::endl;

    // result of the product
    std::vector<data_type> b(A.rows());

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    clock.start();
    b = A * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();
    clock.start();
    b = A * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, row-wise storage: "
              << clock << std::endl;

    // Column wise storage
    algebra::Matrix<data_type, algebra::StorageOrder::COL_WISE> M;
    M.read_mtx(filename);

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    clock.start();
    b = M * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();
    clock.start();
    b = M * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, column-wise storage: "
              << clock << std::endl;





    /////////////////////////////// TEST THE NORM COMPUTATION ///////////////////////////////

    // NORM ONE
    std::cout << "\n\nTESTING NORM ONE COMPUTATION" << std::endl;

    double norm_one;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_one = A.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with uncompressed, row-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_one = A.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with compressed, row-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_one = M.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with uncompressed, column-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_one = M.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with compressed, column-wise storage: " << norm_one << "; " << clock << std::endl;


    // NORM INFINITY
    std::cout << "\n\nTESTING NORM INFINITY COMPUTATION" << std::endl;

    double norm_inf;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_inf = A.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with uncompressed, row-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_inf = A.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with compressed, row-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_inf = M.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with uncompressed, column-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_inf = M.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with compressed, column-wise storage: " << norm_inf << "; " << clock << std::endl;


    // FROBENIUS NORM
    std::cout << "\n\nTESTING FROBENIUS NORM COMPUTATION" << std::endl;

    double norm_fro;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_fro = A.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with uncompressed, row-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_fro = A.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with compressed, row-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_fro = M.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with uncompressed, column-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_fro = M.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Norm infinity with compressed, column-wise storage: " << norm_fro << "; " << clock << std::endl;






/*
    /////////////////////////////// DATA ARE COMPLEX ///////////////////////////////
    std::cout << "COMPLEX DATA\n\n\n"  << std::endl; 

    using data_type = std::complex<double>;

    // reads the required matrix
    algebra::Matrix<data_type> A;
    A.read_mtx(filename);

    // generates the vector to multiply
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    std::vector<data_type> v(A.columns());
    for (std::size_t i=0; i<A.columns(); ++i)
        v[i]=std::complex(uniform(generator),uniform(generator));





    /////////////////////////////// TEST THE PRODUCT MATRIX - VECTOR ///////////////////////////////
    std::cout << "\n\nTESTING MATRIX - VECTOR PRODUCT" << std::endl;

    // result of the product
    std::vector<data_type> b(A.rows());

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    clock.start();
    b = A * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();
    clock.start();
    b = A * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, row-wise storage: "
              << clock << std::endl;

    // Column wise storage
    algebra::Matrix<data_type, algebra::StorageOrder::COL_WISE> M;
    M.read_mtx(filename);

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    clock.start();
    b = M * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();
    clock.start();
    b = M * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, column-wise storage: "
              << clock << std::endl;





    /////////////////////////////// TEST THE NORM COMPUTATION ///////////////////////////////

    // NORM ONE
    std::cout << "\n\nTESTING NORM ONE COMPUTATION" << std::endl;

    double norm_one;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_one = A.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with uncompressed, row-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_one = A.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with compressed, row-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_one = M.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with uncompressed, column-wise storage: " << norm_one << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_one = M.norm<algebra::NormType::One>();
    clock.stop();
    std::cout << "Norm one with compressed, column-wise storage: " << norm_one << "; " << clock << std::endl;


    // NORM INFINITY
    std::cout << "\n\nTESTING NORM INFINITY COMPUTATION" << std::endl;

    double norm_inf;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_inf = A.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with uncompressed, row-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_inf = A.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with compressed, row-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_inf = M.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with uncompressed, column-wise storage: " << norm_inf << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_inf = M.norm<algebra::NormType::Infinity>();
    clock.stop();
    std::cout << "Norm infinity with compressed, column-wise storage: " << norm_inf << "; " << clock << std::endl;


    // FROBENIUS NORM
    std::cout << "\n\nTESTING FROBENIUS NORM COMPUTATION" << std::endl;

    double norm_fro;

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    A.uncompress();

    clock.start();
    norm_fro = A.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with uncompressed, row-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    A.compress();

    clock.start();
    norm_fro = A.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with compressed, row-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    M.uncompress();

    clock.start();
    norm_fro = M.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Frobenius norm with uncompressed, column-wise storage: " << norm_fro << "; " << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    M.compress();

    clock.start();
    norm_fro = M.norm<algebra::NormType::Frobenius>();
    clock.stop();
    std::cout << "Norm infinity with compressed, column-wise storage: " << norm_fro << "; " << clock << std::endl;
*/

    return 0;
} 
