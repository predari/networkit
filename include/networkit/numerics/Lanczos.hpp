/*
 * Lanczos.h
 *
 *  Created on: 30.03.2021
 *      Author: Maria Predari (predarim@hu-berlin.de)
 */


#ifndef NETWORKIT_NUMERICS_LANCZOS_HPP_
#define NETWORKIT_NUMERICS_LANCZOS_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/algebraic/SymTriMatrix.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>
#include <networkit/numerics/Preconditioner/IdentityPreconditioner.hpp>


namespace NetworKit {

/**
 * @ingroup numerics
 * @brief   Implementation of Lanczos Triagonalization.
 * Including algorithms for computing eigenvalues of symmetric diagonal matrices,
 * and retrieving eigenvectors throught inverse power iteration.
 * @param   matrix  input matrix
 * @param   k       number of largest eigenvalues to compute
 * @param   steps   maximum steps for the iteration
 * @tparam  T       element data type
 * @return  list of eigenvalues
 */
template<class Matrix, typename T>
class Lanczos : public Algorithm {
    
public:
    Lanczos(){};
    Lanczos(const Matrix & A);
    Lanczos(const Matrix & A, const int k, int a, int skip, bool double_precision = true, const T epsilon = 1e-5);
    
    ~Lanczos() override = default;

    void setup(const Matrix &A, const int k, int a, int skip, bool double_precision = true, const T epsilon = 1e-5);

    void setup(const Matrix &A);
    
    void run() override;

     std::vector<T> getkEigenvalues();
     T getEigenvalue(int i);

    void computekEigenvectors();
    std::vector<Vector> getkEigenvectors();

    void computeEigenvector(int i);
    Vector getEigenvector(int i);
    bool checkEigenvectors();
    
    template <typename L>
    void forAllEigenvalues(L handle) {
        assureFinished();
        for (typename std::vector<T>::iterator it = eigenvalues.begin(); it != eigenvalues.end(); ++it) {
            handle(*it);
        }
    }
    
    template <typename L>
    void forAllEigenvectors(L handle) {
        assureFinished();
        for (typename std::vector<Vector>::iterator it = eigenvectors.begin(); it != eigenvectors.end(); ++it) {
            handle(*it);
        }
    }
  
    template <typename L>
    void parallelForAllEigenvectors(L handle) {
        assureFinished();
#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(eigenvectors.size()); ++i) {
            handle(eigenvectors[i]);
        }
    }
    
   
  
protected:
    Matrix A;
    int k;
    int a;
    int b;
    int skip;
    T epsilon;
    int steps = 0;
    std::vector<T> eigenvalues;
    std::vector<Vector> eigenvectors;
    
private:
    /**
     * @brief   Performing lanczos iterations.
     * @details creates the symmetric tridiagonal matrix for a required number of steps 
     * @param   number of steps to run the iterations
     * @param   epsilon precision threshold
     * @return  list of eigenvalues
     */
    std::vector<T> lanczos_eigen(int steps);

    /**
     * @brief   Removes spurious eigenvalues
     *
     * @param   number of steps to run the iterations
     * @param   epsilon precision for the check of spurious eigenvalues
     * @return  list of eigenvalues
     */
    std::vector<T> lanczos_no_spurious(SymTriMatrix<T> &tridiag, const T epsilon = 1e-3);

    /**
     * @brief   Calculating eigenvalues for symmetric tridiagonal matrices.
     * @details Reinsch, C. H. (1973). Algorithm 464: Eigenvalues of a Real, Symmetric, 
     * Tridiagonal Matrix. Communications of the ACM, 16(11), 689.
     * 
     * @param   matrix  symmetric tridiagonal matrix to decompose
     * @param   epsilon precision threshold
     * @return  list of eigenvalues
     */
    std::vector<T> tqlrat_eigen(const SymTriMatrix<T> &matrix);
    /**
     * @brief   QR eigendecomposition for symmetric tridiagonal matrices.
     * 
     * @param   matrix  symmetric tridiagonal matrix to decompose
     * @param   epsilon precision threshold
     * @return  list of eigenvalues
     */
    std::vector<T> qr_eigen(const SymTriMatrix<T> &matrix);

    /**
     * @brief   Computing eigenvectors for a given eigenvalue.
     * @details  Implements inverse iteration (power method with a conditioned matrix).
     * If the inverse of A is A', then, 1/eigenvalue is an eigenvalue of A'.
     * Shift the spectrum of A by eigenvalue (now the largest eigenvalue of A').
     * Use the power method on shifted inverse: A'x_i = x_i+1, --> x_i = Ax_i+1 
     * and solve with CG to find x_i+1.
     *
     * @param   eigenvalue 
     * @return  eigenvector
     */
    
    Vector inverse_iteration(T eigenvalue);
    bool inverse_iterationB(T eigenvalue);


};

  template<class Matrix, typename T> Lanczos<Matrix, T> :: Lanczos(const Matrix &A) {
    this->A = A;
    k = 0;
    a = 0;
    skip = 0;
    b = 0;
    epsilon = 0.;
  }


    template<class Matrix, typename T> Lanczos<Matrix, T> :: Lanczos(const Matrix &A, const int k, int a, int skip, bool double_precision, const T epsilon) {

    this->A = A;
    this->k = k;
    int n = A.numberOfColumns();
    if ( (this->k < 1) || (this->k > n) )
      throw std::runtime_error("Error, k = 0 or > n.");
   
    this->a = a;
    this->skip = skip;
    if (double_precision)
      this->b = 64;
    else
      this->b = 8;
    this->epsilon = epsilon;
    
    eigenvalues.resize(k, 0.);
    eigenvectors.resize(k);
    
  }

    template<class Matrix, typename T> void Lanczos<Matrix, T> :: setup(const Matrix &A) {
        
        if ( A == 0 )
            throw std::runtime_error("Error, empty matrix A.");
        this->A = A; 
        int n = A.numberOfColumns();
        if ( (k < 1) || (k > n) )
            throw std::runtime_error("Error, k = 0 or > n.");

}


    
    template<class Matrix, typename T> void Lanczos<Matrix, T> :: setup(const Matrix &A, const int k, int a, int skip, bool double_precision, const T epsilon) {
        
    if ( A == 0 )
      throw std::runtime_error("Error, empty matrix A.");
    this->A = A; 
    this->k = k;
    int n = A.numberOfColumns();
    if ( (this->k < 1) || (this->k > n) )
      throw std::runtime_error("Error, k = 0 or > n.");
    
    this->a = a;
    this->skip = skip;
    if (double_precision)
      this->b = 64;
    else
      this->b = 8;
    this->epsilon = epsilon;
        
    eigenvalues.resize(k, 0.);
    eigenvectors.resize(k);
}
  
  template<class Matrix, typename T> void Lanczos<Matrix, T> :: run () {
  
    std::cout << "*** Running Lanczos ***" << std::endl;
    std::cout << "*** Eigenvalue count k =  " << k <<  std::endl;
    std::cout << "*** Precision epsilon  =  " << epsilon  << std::endl;
    assert(eigenvalues.size() <= k);
    for (int steps = a * k + 1; steps < b * k; steps += skip) {
      eigenvalues = lanczos_eigen(steps);
    }
    std::cout << "*** Lanczos iterations = " <<  this->steps << std::endl;
    std::cout << "*** run status: success ***" << std::endl;

    hasRun = true;
}

  template<class Matrix, typename T> inline std::vector<T> Lanczos<Matrix, T> :: getkEigenvalues() {
    assureFinished();
    if(eigenvalues.size() < k)
        WARN("Lanczos algorithm returns less than k eigenvalues due to spurious behaviour.");
    return eigenvalues;
}

  template<class Matrix, typename T> inline T Lanczos<Matrix, T>::getEigenvalue(int i) {
    assureFinished();
    assert(i < eigenvalues.size());
    return eigenvalues[i];
}

    template<class Matrix, typename T> void Lanczos<Matrix, T>::computekEigenvectors() {
        // recomputation of eigenvectors for every call to computekEigenvectors.
        assureFinished();
        int n = A.numberOfColumns();
        assert(eigenvectors.size() <= k);
        for (int t = 0; t < eigenvectors.size(); ++t) {
            eigenvectors[t] = Vector(inverse_iteration(eigenvalues[t]));
            assert(eigenvectors[t].getDimension() == n);
        }
    }
  

    template<class Matrix, typename T> inline  std::vector<Vector> Lanczos<Matrix, T>::getkEigenvectors() {
        assureFinished();
        assert(eigenvectors.size() <= k);
        // not checking if all eigenvectors are already calculated
        return eigenvectors;
}

    template<class Matrix, typename T> void Lanczos<Matrix, T>::computeEigenvector(int i) {
        assureFinished();
        assert(eigenvectors.size() <= k);
        if (eigenvectors[i].getDimension() == 0)
            eigenvectors[i] = Vector(inverse_iteration(eigenvalues[i]));
    }
  
  
    template<class Matrix, typename T> inline Vector Lanczos<Matrix, T>::getEigenvector(int i) {
        assureFinished();
        assert(eigenvectors.size() <= k);
        assert(i < eigenvectors.size());
        if (eigenvectors[i].getDimension() == 0)
            eigenvectors[i] = Vector(inverse_iteration(eigenvalues[i]));
        return eigenvectors[i];
        
    }

    template<class Matrix, typename T> bool Lanczos<Matrix, T>::checkEigenvectors() {
        assureFinished();
        const T equal_tol  = 1e-4;
        int n = A.numberOfColumns();
        assert(eigenvectors.size() <= k);
        for (int t = 0; t < eigenvectors.size(); ++t) {
            assert(eigenvectors[t].getDimension() == n);
            Vector r1 = A * eigenvectors[t];
            Vector r2 = eigenvalues[t] * eigenvectors[t];
            if (r1.getDimension() != r2.getDimension() || r1.isTransposed() != r2.isTransposed())
                return false;
            for (index i = 0; i < r1.getDimension(); i++) {
                if (!(fabs(r1[i] - r2[i]) < equal_tol))
                    return false;
            }
        }
        return true;
    }


    
  /* TODO: getBasisCSRMatrix() */



  
  
    template <class Matrix, typename T> std::vector<T> Lanczos<Matrix, T> :: lanczos_eigen(int steps) {

        int n = A.numberOfColumns();
        assert(n > 0 && n == A.numberOfRows());
        assert(steps > 2 * k);
  
        SymTriMatrix<T> tridiag(steps);
  
        std::vector< Vector > basis;

        Vector r(n, 0);
        r[0] = 1; // initialize a "random" vector
        T beta = r.length(); // l2 norm
        
        for (int t = 0; t < steps; ++t) {
            if (t > 0) {
                tridiag.beta(t - 1) = beta;
            }
            r = r * (1 / beta);
            basis.push_back(r);
            r = A * r;
            T alpha  = r.transpose() * basis[t];
            r += -alpha * basis[t];
            if (t > 0) {
                r += -beta * basis[t - 1];
            }
            tridiag.alpha(t) = alpha;
            beta = r.length();
        }

        this->steps = steps;  
        // epsilon precision for the spurious is different, default 1e-03
        return lanczos_no_spurious(tridiag);
}

  
    template<class Matrix, typename T>
    std::vector<T> Lanczos<Matrix,T> :: lanczos_no_spurious(SymTriMatrix<T> &tridiag, const T epsilon) {
               
        assert(tridiag.size() > 0);
        std::vector<T> e = tqlrat_eigen(tridiag);
        // for (int i=0; i < e.size(); i++)
        //     std::cout << e[i] << " ";
        // std::cout << std::endl;

        tridiag.remove_forward(0);
        std::vector<T> test_e = tqlrat_eigen(tridiag);
        // for (int i=0; i < test_e.size(); i++)
        //     std::cout << test_e[i] << " ";
        // std::cout << std::endl;

        std::vector<T> eigen;

        int i = 0;
        int j = 0;
        while (j <= e.size()) { // scan through one position beyond the end of the list
            if (j < e.size() && std::abs(e[j] - e[i]) < epsilon) {
                j++;
                continue;
            }
            // simple eigenvalues not in test set are preserved
            // multiple eigenvalues are only preserved once
            if (j - i > 1 || [&] (){
                                 typename std::vector<T>::iterator first = test_e.begin();
                                 typename std::vector<T>::iterator last = test_e.end();
                                 while (first != last) {
                                     if (T(std::abs(*first - e[i])) < epsilon) {
                                         return first;
                                     }
                                     ++first;
                                 }
                                 return last;
                                 
                             }() == test_e.end()) {
                eigen.push_back(e[i]);
            }
            i = j++;
        }
        std::sort(eigen.rbegin(), eigen.rend());

        eigen.resize(std::min((int)eigen.size(), this->k));
        eigenvectors.resize(std::min((int)eigen.size(), this->k));
        return eigen;
    }


    template<class Matrix, typename T>
    std::vector<T> Lanczos<Matrix,T> :: tqlrat_eigen(const SymTriMatrix<T> &matrix) {
        int n = matrix.size();
        std::vector<T> d(matrix.alpha_data(), matrix.alpha_data() + n);
        std::vector<T> e2(n, 0);
        for (int i = 0; i < n - 1; ++i) {
            e2[i] = matrix.beta(i) * matrix.beta(i);
        }
        T b(0), b2(0), f(0);
        for (int k = 0; k < n; ++k) {
            T h = epsilon * epsilon * (d[k] * d[k] + e2[k]);
            if (b2 < h) {
                b = sqrt(h);
                b2 = h;
            }
            int m = k;
            while (m < n && e2[m] > b2) {
                ++m;
            }
            if (m == n) {
                --m;
        }
            if (m > k) {
                do {
                    T g = d[k];
                    T p2 = sqrt(e2[k]);
                    h = (d[k + 1] - g) / (2.0 * p2);
                    T r2 = sqrt(h * h + 1.0);
                    d[k] = h = p2 / (h < 0.0 ? h - r2 : h + r2);
                    h = g - h;
                    f = f + h;
                    for (int i = k + 1; i < n; ++i) {
                        d[i] -= h;
                    }
                    h = g = std::abs(d[m] - 0.0) < epsilon ? b : d[m];
                    T s2 = 0.0;
                    for (int i = m - 1; i >= k; --i) {
                        p2 = g * h;
                        r2 = p2 + e2[i];
                        e2[i + 1] = s2 * r2;
                        s2 = e2[i] / r2;
                        d[i + 1] = h + s2 * (h + d[i]);
                        g = d[i] - e2[i] / g;
                        if (std::abs(g - 0.0) < epsilon) {
                            g = b;
                        }
                        h = g * p2 / r2;
                    }
                    e2[k] = s2 * g * h;
                    d[k] = h;
                } while (e2[k] > b2);
            }
            h = d[k] + f;
            int j;
            for (j = k; j > 0; --j) {
                if (h >= d[j - 1]) {
                    break;
                }
                d[j] = d[j - 1];
            }
            d[j] = h;
        }

        return d;
    }


    template<class Matrix, typename T>
    std::vector<T> Lanczos<Matrix,T> :: qr_eigen(const SymTriMatrix<T> &matrix) {


        SymTriMatrix<T> tridiag = matrix;
        int n = tridiag.size();

        tridiag.resize(n + 1);
        tridiag.alpha(n) = 0;
        tridiag.beta(n - 1) = 0;
        for (int i = 0; i < n - 1; ++i) {
            tridiag.beta(i) = tridiag.beta(i) * tridiag.beta(i);
        }
        bool converged = false;
        while (!converged) {
            T diff(0);
            T u(0);
            T ss2(0), s2(0); // previous and current value of s^2
            for (int i = 0; i < n; ++i) {
                T gamma = tridiag.alpha(i) - u;
                T p2 = T(std::abs(1 - s2)) < epsilon ? (1 - ss2) * tridiag.beta(i - 1) : gamma * gamma / (1 - s2);
                if (i > 0) {
                    tridiag.beta(i - 1) = s2 * (p2 + tridiag.beta(i));
                }
                ss2 = s2;
                s2 = tridiag.beta(i) / (p2 + tridiag.beta(i));
                u = s2 * (gamma + tridiag.alpha(i + 1));
                // update alpha
                T old = tridiag.alpha(i);
                tridiag.alpha(i) = gamma + u;
                diff = std::max(diff, T(std::abs(old - tridiag.alpha(i))));
            }
            if (diff < epsilon) {
                converged = true;
            }
        }
        
        return std::vector<T>(tridiag.alpha_data(), tridiag.alpha_data() + n);
}


  template<class Matrix, typename T>
  Vector Lanczos<Matrix,T> :: inverse_iteration(T eigenvalue) {


    int n = A.numberOfColumns();
    Vector eigenvector = Vector(n);
    for (index i = 0; i < n; ++i) {
        eigenvector[i] = 2.0 * Aux::Random::real() - 1.0;
    }

    const double solver_tol = 1e-5;
    Vector old;
    count numIterations = 0;
    const auto lambdaI = Matrix::diagonalMatrix(eigenvalue*Vector(n,1.0));
    Matrix M = A - lambdaI;
    
    do {
        old = eigenvector;
        ConjugateGradient<Matrix, IdentityPreconditioner<Matrix>> cg(solver_tol);
        // TODO: must solve for tridiagonal matrix instead!!! Use basis from Lanczos to compute eigenvectors of A.
        cg.setup(M);
        cg.solve(old, eigenvector);
        eigenvector /= eigenvector.length();
      
        numIterations++;
    } while ((eigenvector - old).length() > 1e-6 && numIterations < 1500);
    
    return eigenvector;
}



  template<class Matrix, typename T>
  bool Lanczos<Matrix,T> :: inverse_iterationB(T eigenvalue) {

    double halt = 1e-4;
    int n = A.numberOfColumns();
    int rows = A.numberOfColumns();
    int iterations = rows;

    const auto lambdaI = Matrix::diagonalMatrix(eigenvalue*Vector(n,1.0));
    Matrix _A = A - lambdaI;
    Matrix scratch_A;
    Vector d;
    double r_inf;
    Vector x = Vector(n);

    Vector u = Vector(n);
    for (index i = 0; i < n; ++i) {
      u[i] = 2.0 * Aux::Random::real() - 1.0;
    }
    
    //u.set_WiP ();
    // _Ax_n = x_n+1 -> Ax_n+1=x_n, factor with QR
    while (true) {
        scratch_A = _A;
        //scratch_A.copy ();
        ConjugateGradient<Matrix, IdentityPreconditioner<Matrix>> cg(1e-5);
        cg.setup(scratch_A);
        cg.solve(u, x);
        x /= x.length();
        u = x;
        d = _A * u;


        r_inf = DBL_MIN;
        for (int i = 0; i < rows; ++i)
            if (r_inf < fabs (d[i]))
                r_inf = fabs (d[i]);
        
#if 0

        Vector Rayleigh = A * u;
       	double quotient = u.transpose () * Rayleigh;;
        double r = ((A * u - eigenvalue * u).length());
        double R = eigenvalue - quotient;
#endif

        if (r_inf <= halt)
            break;
        
        --iterations;
        if (iterations < 0)
                return false;
    }
    
    return true;
}

  

  

  
} /* namespace NetworKit */
  
#endif // NETWORKIT_NUMERICS_LANCZOS_HPP_
