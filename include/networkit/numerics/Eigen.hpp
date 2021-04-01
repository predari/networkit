/*
 * Eigen.h
 *
 *  Created on: 30.03.2021
 *      Author: Maria Predari (predarim@hu-berlin.de)
 */

#ifndef NETWORKIT_NUMERICS_EIGEN_HPP_
#define NETWORKIT_NUMERICS_EIGEN_HPP_

#include <networkit/algebraic/Vector.hpp>
#include <networkit/auxiliary/Timer.hpp>
// #include <networkit/algebraic/CSRMatrix.hpp>
namespace NetworKit {

/**
 * @brief A symmetric tridiagonal matrix.
 * @details
 * 
 * @tparam T Element data type.
 */
template <typename T>
class symm_tridiag_matrix {
public:
    symm_tridiag_matrix(int n)
        : alpha_(n), beta_(n - 1) {}

    int size() const { return alpha_.size(); }
    void resize(int n);
    void remove_forward(int i);
    void remove_backward(int i);

    const T &alpha(int i) const { return alpha_[i]; }
    const T &beta(int i) const { return beta_[i]; }

    T &alpha(int i) { return alpha_[i]; }
    T &beta(int i) { return beta_[i]; }

    const T *alpha_data() const { return alpha_.data(); }
    const T *beta_data() const { return beta_.data(); }

    T *alpha_data() { return alpha_.data(); }
    T *beta_data() { return beta_.data(); }

private:
  std::vector<T> alpha_; /**< main diagonal entries */
  std::vector<T> beta_; /**< diagonal entries below or above the main diagonal */
};

template <typename T>
void symm_tridiag_matrix<T>::resize(int n) {
    alpha_.resize(n);
    beta_.resize(n - 1);
}

template <typename T>
void symm_tridiag_matrix<T>::remove_forward(int i) {
    assert(i < size() - 1);
    alpha_.erase(alpha_.begin() + i);
    beta_.erase(beta_.begin() + i);
}

template <typename T>
void symm_tridiag_matrix<T>::remove_backward(int i) {
    assert(i > 0);
    alpha_.erase(alpha_.begin() + i);
    beta_.erase(beta_.begin() - 1 + i);
}


/**
 * @brief   Lanczos algorithm for eigendecomposition.
 * 
 * @param   matrix  CSR matrix to decompose
 * @param   k       number of largest eigenvalues to compute
 * @param   steps   maximum steps for the iteration
 * @tparam  T       matrix element data type
 * @return  list of eigenvalues
 */


  
  template <class Matrix, typename T = double>
  std::vector<T> lanczos_eigen(const Matrix& matrix, int k, int steps) {


  int n = matrix.numberOfColumns();
  assert(n > 0 && n == matrix.numberOfRows());
  assert(steps > 2 * k);
  Aux::Timer timer;
  
  symm_tridiag_matrix<T> tridiag(steps);
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
    r = matrix * r;
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

  return lanczos_no_spurious(tridiag, k);
  
  //std::vector<T> v = lanczos_no_spurious(tridiag, k);
  // Vector f(v);
  // f.forElements([&](const double u) {
  // 		  std::cout << "node = " << u   << "\n";
  // 		});

  // for ( int i = 0; i< basis.size(); i++ ) {
  //   std::cout << " eigenvector #" << i <<  "[ ";
  //   basis[i].forElements([&](const double u) {
  // 			   std::cout  << u << " ";
  // 		});
  //     std::cout << " ]\n";
  // }
  // return f;
  ////return Vector(lanczos_no_spurious(tridiag, k));
  
}

template <typename T>
std::vector<T> lanczos_no_spurious(symm_tridiag_matrix<T> &tridiag, int k, const T epsilon = 1e-3) {
    assert(tridiag.size() > 0);
    Aux::Timer timer;
    timer.start(); 

    std::vector<T> eigen = tqlrat_eigen(tridiag);
    tridiag.remove_forward(0);
    std::vector<T> test_eigen = tqlrat_eigen(tridiag);
    std::vector<T> result;

    int i = 0;
    int j = 0;
    while (j <= eigen.size()) { // scan through one position beyond the end of the list
        if (j < eigen.size() && std::abs(eigen[j] - eigen[i]) < epsilon) {
            j++;
            continue;
        }
        // simple eigenvalues not in test set are preserved
        // multiple eigenvalues are only preserved once
        if (j - i > 1 || [&] (){
			   std::vector<double>::iterator first = test_eigen.begin();
			   std::vector<double>::iterator last = test_eigen.end();
			   while (first != last) {
			     if (T(std::abs(*first - eigen[i])) < epsilon) {
			       return first;
			     }
			     ++first;
			   }
			   return last;
			   
			 }() == test_eigen.end()) {
            result.push_back(eigen[i]);
        }
        i = j++;
    }
    std::sort(result.rbegin(), result.rend());
    result.resize(std::min((int)result.size(), k));

    timer.stop(); 
    std::cout << "spurious removal time: " << timer.elapsedMicroseconds() / 1e6 << " sec" << std::endl;
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
template <typename T>
std::vector<T> tqlrat_eigen(const symm_tridiag_matrix<T> &matrix, const T epsilon = 1e-8) {

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
template <typename T>
std::vector<T> qr_eigen(const symm_tridiag_matrix<T> &matrix, const T epsilon = 1e-8) {

  Aux::Timer timer;
  timer.start(); 
    symm_tridiag_matrix<T> tridiag = matrix;
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

  
} /* namespace NetworKit */
  
#endif // NETWORKIT_NUMERICS_EIGEN_HPP_
