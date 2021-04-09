/*
 * Lanczos.h
 *
 *  Created on: 30.03.2021
 *      Author: Maria Predari (predarim@hu-berlin.de)
 */

#include <iostream>
// fabs
#include <cmath>


#ifndef NETWORKIT_NUMERICS_EIGEN_HPP_
#define NETWORKIT_NUMERICS_EIGEN_HPP_

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
 * @brief   Implementation of Lanczos.
 * @param   matrix  CSR matrix to decompose
 * @param   k       number of largest eigenvalues to compute
 * @param   steps   maximum steps for the iteration
 * @tparam  T       matrix element data type
 * @return  list of eigenvalues
 */

  
 
template<class Matrix, typename T>
class Lanczos : public Algorithm {
    
public:
    Lanczos(){};
    Lanczos(const Matrix & A);
    Lanczos(const Matrix & A, const int k, int a, int skip, bool double_precision = false);
    
    ~Lanczos() override = default;

    void setup(const int k, int a, int skip, bool double_precision = false);

    
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
        for (typename std::vector<T>::iterator it = eigen.begin(); it != eigen.end(); ++it) {
            handle(*it);
        }
    }
    
    template <typename L>
    void forAllEigenvectors(L handle) {
        assureFinished();
        for (typename std::vector<Vector>::iterator it = basis.begin(); it != basis.end(); ++it) {
            handle(*it);
        }
    }
  
    template <typename L>
    void parallelForAllEigenvectors(L handle) {
        assureFinished();
#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(basis.size()); ++i) {
            handle(basis[i]);
        }
    }
    
   
  
protected:
    Matrix A;
    int k;
    int a;
    int b;
    int skip;
    T tolerance;
    std::vector<T> eigen;
    std::vector<Vector> basis;
    
private:
    std::vector<T> lanczos_eigen(int steps);
    std::vector<T> lanczos_no_spurious(SymTriMatrix<T> &tridiag, const T epsilon = 1e-3);
  
    /* algorithm to compute the eigenvalues of a symmetric tridiagonal matrix */ 
    std::vector<T> tqlrat_eigen(const SymTriMatrix<T> &matrix, const T epsilon = 1e-8);
    /* algorithm to compute the eigenvalues of a symmetric tridiagonal matrix (qr-based)*/ 
    std::vector<T> qr_eigen(const SymTriMatrix<T> &matrix, const T epsilon = 1e-8);
    /* inverse power method to compute eigenvector */
    Vector find_eigen_vector(T eigenvalue);
    bool find_eigen_vector2 (T eigenvalue);
  

};

  template<class Matrix, typename T> Lanczos<Matrix, T> :: Lanczos(const Matrix &A) {
    this->A = A;
    this->k = 0;
    this->a = 0;
    this->skip = 0;
    this->b = 0;
  }


  template<class Matrix, typename T> Lanczos<Matrix, T> :: Lanczos(const Matrix &A, const int k, int a, int skip, bool double_precision) {

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
    
    eigen.resize(k, 0.);
    basis.resize(k);
    
  }
  
  template<class Matrix, typename T> void Lanczos<Matrix, T> :: setup(const int k, int a, int skip, bool double_precision) {

    if ( A == 0 )
      throw std::runtime_error("Error, empty matrix A.");
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
    
    eigen.resize(k, 0.);
    basis.resize(k);

    //distanceToTarget.assign(G.upperNodeIdBound(), none);
    //distanceFromSource.assign(G.upperNodeIdBound(), none);


}
  
  template<class Matrix, typename T> void Lanczos<Matrix, T> :: run () {
  
    std::cout << "*** running CPU Lanczos ***" << std::endl;
    std::cout << "*** k =  " << k << " a = " << a << " skip = " << skip << " b = " <<  b  << std::endl;
    
    for (int steps = a * k + 1; steps < b * k; steps += skip) {
      eigen = lanczos_eigen (steps);
    }  
    std::cout << "*** run status: success ***" << std::endl;
    hasRun = true;
}

  template<class Matrix, typename T> inline std::vector<T> Lanczos<Matrix, T> :: getkEigenvalues() {
    assureFinished();
    return eigen;
}

  template<class Matrix, typename T> inline T Lanczos<Matrix, T>::getEigenvalue(int i) {
    assureFinished();
    assert(i < eigen.size());
    return eigen[i];
}

  template<class Matrix, typename T> inline  void Lanczos<Matrix, T>::computekEigenvectors() {
    assureFinished();
    int n = A.numberOfColumns();
    std::cout << " n = " << n << " basis size = " << basis.size() << std::endl;
    assert(basis.size() > 0);
    for (int t = 0; t < basis.size(); ++t) {
      if (basis[t].getDimension() <= 0) {
	Vector r = find_eigen_vector(eigen[t]);
	//find_eigen_vector2(eigen[t]);
	std::cout << "*** computing eigenvector " << t << " based on eigenvalue " << eigen[t] <<  " : success ***" << std::endl;
	assert(r.getDimension() == n);
	basis[t] = Vector(r);
      }
    }
  }

  template<class Matrix, typename T> inline bool Lanczos<Matrix, T>::checkEigenvectors() {
    assureFinished();
    int n = A.numberOfColumns();
    std::cout << " n = " << n << " basis size = " << basis.size() << std::endl;
    assert(basis.size() > 0);
    for (int t = 0; t < basis.size(); ++t) {
      assert(basis[t].getDimension() == n);
      Vector r1 = A*basis[t];
      Vector r2 = eigen[t] * basis[t];
      if ( r1 == r2)
	std::cout << "*** Eigenvectors are computed successfully! ***" << std::endl;
      else {
	std::cout << "*** WARNING!!!! Eigenvector " << t << " is wrong! ***" << std::endl;
	return false;
      }
    }
    return true;
  }

  

  template<class Matrix, typename T> inline  std::vector<Vector> Lanczos<Matrix, T>::getkEigenvectors() {
    assureFinished();
    // if(basis.size() > 0 && basis[0].getDimension() == n) return basis;
    // else computekEigenvectors();
    return basis;
}

  template<class Matrix, typename T> inline  void Lanczos<Matrix, T>::computeEigenvector(int i) {
    assureFinished();
    assert(basis.size() > 0);
    if (basis[i].getDimension() == 0)
      basis[i] = Vector(find_eigen_vector(eigen[i]));
    return;
}



  
  
  template<class Matrix, typename T> inline Vector Lanczos<Matrix, T>::getEigenvector(int i) {
    assureFinished();
    assert(basis.size() > 0);
    assert(i < basis.size());
    if (basis[i].getDimension() == 0)
      basis[i] = Vector(find_eigen_vector(eigen[i]));
    //return basis[i];
    return Vector(find_eigen_vector(eigen[i]));
    
}

  /* TODO: getBasisCSRMatrix() */



  
  
  template <class Matrix, typename T>
  std::vector<T> Lanczos<Matrix, T> :: lanczos_eigen(int steps) {

  int n = A.numberOfColumns();
  assert(n > 0 && n == A.numberOfRows());
  assert(steps > 2 * k);
  Aux::Timer timer;
  SymTriMatrix<T> tridiag(steps);
  std::cout << " ** 1st : size of tridiagonal : " << tridiag.size() << std::endl;
  
  std::vector< Vector > basis;

  Vector r(n, 0);
  r[0] = 1; // initialize a "random" vector
  T beta = r.length(); // l2 norm
  timer.start();

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
  timer.stop();
  std::cout << "CPU Lanczos iterations: " << steps << std::endl;
  double timeElapsed = timer.elapsedMicroseconds() / 1e6;
  std::cout << "CPU Lanczos time: " << timeElapsed << " sec" << std::endl;

  std::cout << " ** 2nd : size of tridiagonal : " << tridiag.size() << std::endl;
  std::cout << " ** 3rd : size of basis : " << basis.size() << " of size vectors : " << basis[0].length() << std::endl;
  std::cout << " Calling lanczos_no_spurious "  << std::endl;

  std::vector<T> result = lanczos_no_spurious(tridiag, k);
  return result;
  //return lanczos_no_spurious(tridiag, k);
  
  
}

  
  template<class Matrix, typename T>
  std::vector<T> Lanczos<Matrix,T> :: lanczos_no_spurious(SymTriMatrix<T> &tridiag, const T epsilon) {
    assert(tridiag.size() > 0);
    // Aux::Timer timer;
    // timer.start(); 
    std::cout << " ** 4rth : size of tridiagonal : " << tridiag.size() << std::endl;
    std::vector<T> e = tqlrat_eigen(tridiag);
    std::cout << " ** 5th : size of tridiagonal after tqlrat : " << tridiag.size() << std::endl;
    tridiag.remove_forward(0);
    std::cout << " ** 6th : size of tridiagonal after remove_forward : " << tridiag.size() << std::endl;
    std::vector<T> test_e = tqlrat_eigen(tridiag);
    std::cout << " ** 7th : size of tridiagonal after second tqlrat : " << tridiag.size() << std::endl;
    std::vector<T> result;

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
            result.push_back(e[i]);
        }
        i = j++;
    }
    std::sort(result.rbegin(), result.rend());
    result.resize(std::min((int)result.size(), this->k));

    // timer.stop(); 
    // std::cout << "spurious removal time: " << timer.elapsedMicroseconds() / 1e6 << " sec" << std::endl;
    return result;
}

/**
 * @brief   Calculating eigenvalues for symmetric tridiagonal matrices.
 * @details Reinsch, C. H. (1973). Algorithm 464: Eigenvalues of a Real, Symmetric, Tridiagonal Matrix.
 *          Communications of the ACM, 16(11), 689.
 * 
 * @param   matrix  symmetric tridiagonal matrix to decompose
 * @param   epsilon precision threshold
 * @tparam  T       matrix element data type
 * @return  list of eigenvalues
 */
  template<class Matrix, typename T>
  std::vector<T> Lanczos<Matrix,T> :: tqlrat_eigen(const SymTriMatrix<T> &matrix, const T epsilon) {

  Aux::Timer timer;
  timer.start(); 

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


    timer.stop(); 
    std::cout << "TQLRAT time: " << timer.elapsedMicroseconds() / 1e6 << " sec" << std::endl;
    return d;
}

/**
 * @brief   QR eigendecomposition for symmetric tridiagonal matrices.
 * 
 * @param   matrix  symmetric tridiagonal matrix to decompose
 * @param   epsilon precision threshold
 * @tparam  T       matrix element data type
 * @return  list of eigenvalues
 */
  template<class Matrix, typename T>
  std::vector<T> Lanczos<Matrix,T> :: qr_eigen(const SymTriMatrix<T> &matrix, const T epsilon) {

  Aux::Timer timer;
  timer.start(); 
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

    timer.stop();
    std::cout << "QR decomposition time: " << timer.elapsedMicroseconds() / 1e6 << " sec" << std::endl;
    return std::vector<T>(tridiag.alpha_data(), tridiag.alpha_data() + n);
}




  /*
 * Routine to compute eigenvectors for real eigenvalues.  Implements 
 * inverse iteration (power method with a conditioned matrix).
 *
 * If the inverse of A is A', then,
 * 1. 1/µ is an eigenvalue of A'.
 * 2. Shift the spectrum of A by µ (now the largest eigenvalue of A')
 * 3. By the power method: A'x_i = x_i+1, --> x_i = Ax_i+1:
 *	  factor with QR to find x_i+1
 *
 */

  template<class Matrix, typename T>
  Vector Lanczos<Matrix,T> :: find_eigen_vector (T eigenvalue) {

    std::cout << "Finding eigenvector for eigenvalue : " << eigenvalue << std::endl;
    int n = A.numberOfColumns();
    std::cout << "n = " << n << std::endl;
    Vector eigenvector = Vector(n);
    for (index i = 0; i < n; ++i) {
      eigenvector[i] = 2.0 * Aux::Random::real() - 1.0;
    }

    std::cout << "Before iter: Eigenvector : ["<< std::endl;
    for (int i =0; i< eigenvector.getDimension(); i++) {  
      std::cout << eigenvector[i]<< " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Before iter: dimension = " << eigenvector.getDimension() << std::endl;
    
    // const double tol // default 1e-5
    Vector old;
    count numIterations = 0;
    const auto lambdaI = Matrix::diagonalMatrix(eigenvalue*Vector(n,1.0));
    Matrix M = A - lambdaI;
    
    do {
      old = eigenvector;
      // inverse power iteration instead of power iteration: eigenvector = A*old;
      ConjugateGradient<Matrix, IdentityPreconditioner<Matrix>> cg(1e-5);
      // TODO: solving for A!! More performant to sovle for T!!!
      //cg.setup(A);
      cg.setup(M);
      cg.solve(old, eigenvector);
      //
      eigenvector /= eigenvector.length();
      
        numIterations++;
    } while ((eigenvector - old).length() > 1e-6 && numIterations < 10000);//1500);

    std::cout << "Eigenvector: ["<< std::endl;
    for (int i =0; i< eigenvector.getDimension(); i++) {  
      std::cout << eigenvector[i]<< " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "After iter: dimension = " << eigenvector.getDimension() << std::endl;
    return eigenvector;
}



  template<class Matrix, typename T>
  bool Lanczos<Matrix,T> :: find_eigen_vector2 (T eigenvalue) {

    double halt = 1e-4;
    std::cout << "Finding eigenvector for eigenvalue : " << eigenvalue << std::endl;
    int n = A.numberOfColumns();
    std::cout << "n = " << n << std::endl;

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
    std::cout << "PRINT 1 " << eigenvalue << std::endl;
    
    //u.set_WiP ();

    // _Ax_n = x_n+1 -> Ax_n+1=x_n, factor with QR
    while (true) {
      std::cout << "PRINT 2 " << eigenvalue << std::endl;
        scratch_A = _A;
        //scratch_A.copy ();
	ConjugateGradient<Matrix, IdentityPreconditioner<Matrix>> cg(1e-5);
	cg.setup(scratch_A);
	cg.solve(u, x);
	x /= x.length();
	std::cout << "PRINT 3 " << eigenvalue << std::endl;
	u = x;
	assert(_A.numberOfColumns() == _A.numberOfRows());
	assert(_A.numberOfColumns() == u.getDimension());
	d = _A * u;

	r_inf = 1e-5;
	// r_inf = DBL_MIN;
	for (int i = 0; i < rows; ++i)
	  if (r_inf < fabs (d[i]))
	    r_inf = fabs (d[i]);
	
#if 0
	assert(A.numberOfColumns() == A.numberOfRows());
	assert(A.numberOfColumns() == u.getDimension());
        Vector Rayleigh = A * u;
       	double quotient = u.transpose () * Rayleigh;;
        double r = ((A * u - eigenvalue * u).length());
	double R = eigenvalue - quotient;

// 	printf ("Halt conditions %e\t%e\t%e\n", r, halt, A.norm_inf ());
// 	printf ("Rayleigh %e %f\n", quotient / lambda, quotient);
#endif

        if (r_inf <= halt)
	  break;

	--iterations;
	if (iterations < 0)
	  {
	    std::cout << " Didn't converge! " << std::endl;
	    std::cout << " vector u " << std::endl;
	    u.forElements([&](const double &value) {
			    std::cout << value << std::endl;
			  });
	    
	    std::cout << " vector d " << std::endl;
	    d.forElements([&](const double &value) {
			    std::cout << value << std::endl;
			  });
	    
	    return false;
	  }
	  

    }
    std::cout << " vector u " << std::endl;
    u.forElements([&](const double &value) {
		    std::cout << value << std::endl;
		  });

    // std::cout << " vector d " << std::endl;
    // d.forElements([&](const double &value) {
    // 		    std::cout << value << std::endl;
    // 		  });

    
    return true;
}

  

  

  
} /* namespace NetworKit */
  
#endif // NETWORKIT_NUMERICS_EIGEN_HPP_
