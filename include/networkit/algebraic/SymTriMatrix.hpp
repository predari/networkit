/*
 * SymTriMatrix.h
 *
 *  Created on: 30.03.2021
 *      Author: Maria Predari (predarim@hu-berlin.de)
 */

#ifndef NETWORKIT_ALGEBRAIC_SYM_TRI_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_SYM_TRI_MATRIX_HPP_


namespace NetworKit {

/**
 * @brief A symmetric tridiagonal matrix.
 * @details
 * 
 * @tparam T Element data type.
 */
template <typename T>
class SymTriMatrix {
public:
    SymTriMatrix(int n)
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
void SymTriMatrix<T>::resize(int n) {
    alpha_.resize(n);
    beta_.resize(n - 1);
}

template <typename T>
void SymTriMatrix<T>::remove_forward(int i) {
    assert(i < size() - 1);
    alpha_.erase(alpha_.begin() + i);
    beta_.erase(beta_.begin() + i);
}

template <typename T>
void SymTriMatrix<T>::remove_backward(int i) {
    assert(i > 0);
    alpha_.erase(alpha_.begin() + i);
    beta_.erase(beta_.begin() - 1 + i);
}

}

#endif // NETWORKIT_ALGEBRAIC_SYM_TRI_MATRIX_HPP_
