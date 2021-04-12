/*
 * DynBetweenness.h
 *
 *  Created on: 12.08.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#ifndef NETWORKIT_NUMERICS_DYN_LANCZOS_HPP_
#define NETWORKIT_NUMERICS_DYN_LANCZOS_HPP_


#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/numerics/Lanczos.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>


namespace NetworKit {


/**
 * @ingroup numerics
 * Dynamic eigen-updates based on lanczos.
 */

template<class Matrix, typename T>
class DynLanczos : public Lanczos<Matrix, T>, public DynAlgorithm {

  
public:
    /**
     * Creates the object for @a A.
     *
     * @param A The matrix.
     */
    DynLanczos(const Matrix & A, const int k, int a, int skip, bool double_precision = true, const T epsilon = 1e-05);

    ~DynLanczos() override = default;

    /**
     * Runs static eigenpair computation via lanczos algorithm on the initial graph.
     */
    void run() override;
    
    /**
     * Updates the eigenpair computation after an edge event on the graph.
     *
     * @param e The edge event.
     */
    void update(GraphEvent singleEvent) override {
        std::vector<GraphEvent> events{singleEvent};
        updateBatch(events);
    }
    /**
     * Updates the eigenpair computation after a batch of edge events on the graph.
     * Notice: it works only with edge events.
     *
     * @param batch The batch of edge events.
     */
    
    void updateBatch(const std::vector<GraphEvent>& batch) override;

   

    std::vector<T> getUpdatedEigenvalues() ;
    
    T getUpdatedEigenvalue(int i) ;
    
    std::vector<Vector> getUpdatedEigenvectors() ;
  
    Vector getUpdatedEigenvector(int i) ;

    void loadEigenvectors(std::vector<Vector>);
    

private:
  
  Matrix Delta;
  std::vector<T> eigen;
  std::vector<T> eigenOld;
  std::vector<Vector> basis;
  std::vector<Vector> basisOld;
  count affectedEdges = 0;
  double timeDep = 0;
  
};
  
  template<class Matrix, typename T>
  DynLanczos<Matrix, T>::DynLanczos(const Matrix & A, const int k, int a, int skip, bool double_precision, const T epsilon) :
      Lanczos<Matrix,T>(A, k, a, skip, double_precision, epsilon), Delta(Matrix(A.numberOfRows(), A.numberOfColumns())), eigen(k, 0.), basis(k, Vector(A.numberOfRows(), 0.)){ }
   

  
template<class Matrix, typename T> void DynLanczos<Matrix, T> :: run() {

    Lanczos<Matrix,T>::run();
    Lanczos<Matrix,T>::assureFinished();
    Lanczos<Matrix,T>::computekEigenvectors();
    
   
    eigenOld = Lanczos<Matrix,T>::eigenvalues;
    basisOld = Lanczos<Matrix,T>::eigenvectors;
    int n = Lanczos<Matrix,T>::A.numberOfColumns();
    int k = Lanczos<Matrix,T>::k;
    assert(eigenOld.size() <= k);
    assert(basisOld.size() <= k);


    DEBUG("************* after static run: ****************** ");
    std::cout << "[ ";
    for (int i = 0; i < eigenOld.size(); i++)
        std::cout << eigenOld[i] << " ";
    std::cout << "]\n";
    for (int i =0; i< basisOld.size(); i++) {
        assert(basisOld[i].getDimension() == n);
        std::cout << "  ******** Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basisOld[i][j] << " ";
        }
        std::cout << "]\n";
    }

    std::cout << " eigen [ ";
    for (int i = 0; i < eigen.size(); i++)
        std::cout << eigen[i] << " ";
    std::cout << "]\n";


    this->hasRun = true;
}


 
    template<class Matrix, typename T> void DynLanczos<Matrix, T> :: updateBatch(const std::vector<GraphEvent>& batch) {

        INFO("Entering update");
        
        const int n = Lanczos<Matrix,T>::A.numberOfRows();
        assert(!Delta.nnz() && (Delta.numberOfRows() == n));
        
    
        for(auto event : batch){
            node u = event.u;
            node v = event.v;
            if (!(event.type==GraphEvent::EDGE_ADDITION ||
                  event.type==GraphEvent::EDGE_WEIGHT_INCREMENT ||
                  event.type==GraphEvent::EDGE_REMOVAL)) {
                throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
            }
            
            Delta.setValue(u, v, -1.0);
            Delta.setValue(v, u, -1.0);
            Delta.setValue(v, v, Delta(v, v) + event.w);
            Delta.setValue(u, u, Delta(u, u) + event.w);
        }
        

    std::cout << " old[ ";
    for (int i = 0; i < eigenOld.size(); i++)
        std::cout << eigenOld[i] << " ";
    std::cout << "]\n";
    DEBUG(" ******************************** ");
    for (int i =0; i< basisOld.size(); i++) {
        assert(basisOld[i].getDimension() == n);
        std::cout << " *** Old Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basisOld[i][j] << " ";
        }
        std::cout << "]\n";
    }    

        
    
    for(int j = 0; j < Lanczos<Matrix,T>::k; j++) {
      Vector deltaEigenvector(n,0.);
      Vector v(n,0.);
      Vector b = basisOld[j];
      assert(b.getDimension() == n);
      for (int k = 0; k < n; k++) {
          std::cout << b[k] << " ";
      }
         std::cout << "\n";
      for(int i = 0; i < Lanczos<Matrix,T>::k; i++) {
          if (i == j) continue;
          
          
          
          v = Delta * b;
          double alpha = basisOld[i].transpose() * v;
          alpha /=  (eigenOld[j] - eigenOld[i]);
          INFO("eigenOld[j] - eigenOld[i] = ", eigenOld[j] - eigenOld[i]);
          INFO(" delta eigenvector");
          deltaEigenvector = deltaEigenvector + alpha * basisOld[i];
          for (int i= 0; i < n; i++)
              INFO(deltaEigenvector[i]);
          
      }
      double deltaEigenvalue = basisOld[j].transpose() * v;
      INFO("deltaEigenvalue = ", deltaEigenvalue);
      eigen[j] = eigenOld[j] + deltaEigenvalue; 
      basis[j] = basisOld[j] + deltaEigenvector;
    }


    // DEBUG("******************************** ");
    // DEBUG("********* update() ************* ");
    // std::cout << " old[ ";
    // for (int i = 0; i < eigenOld.size(); i++)
    //     std::cout << eigenOld[i] << " ";
    // std::cout << "]\n";
    // DEBUG(" ******************************** ");
    // for (int i =0; i< basisOld.size(); i++) {
    //     assert(basisOld[i].getDimension() == n);
    //     std::cout << " *** Old Eigenvector # " <<  i << "[ ";
    //     for (int j = 0; j < n; j++) {
    //         std::cout << basisOld[i][j] << " ";
    //     }
    //     std::cout << "]\n";
    // }    
    // DEBUG("[DYN] ************************************** ");

    // update for next run
    eigenOld = eigen;
    basisOld = basis;
    Lanczos<Matrix,T>::setup(Lanczos<Matrix,T>::A + Delta);
    Delta = Matrix(Delta.numberOfRows());



    DEBUG(" ********* after update() ************* ");
    std::cout << "new eigen [ ";
    for (int i = 0; i < eigen.size(); i++)
        std::cout << eigen[i] << " ";
    std::cout << "]\n";
    DEBUG(" ******************************** ");
    for (int i =0; i< basis.size(); i++) {
        assert(basis[i].getDimension() == n);
        std::cout << "  *** New Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basis[i][j] << " ";
        }
        std::cout << "]\n";
    }    



    
  }

  template<class Matrix, typename T> inline std::vector<T> DynLanczos<Matrix, T> :: getUpdatedEigenvalues() {
      this->assureFinished();
      if(eigen.size() < Lanczos<Matrix,T>::k)
          WARN("Lanczos algorithm returns less than k eigenvalues due to spurious behaviour.");
      return eigen;
}

  template<class Matrix, typename T> inline T DynLanczos<Matrix, T>::getUpdatedEigenvalue(int i) {
    this->assureFinished();
    assert(i < eigen.size());
    return eigen[i];
}

    template<class Matrix, typename T> inline  std::vector<Vector> DynLanczos<Matrix, T>::getUpdatedEigenvectors() {
        this->assureFinished();
        int k = Lanczos<Matrix,T>::k;
        assert(basis.size() <= k);
        return basis;
    }

    
    template<class Matrix, typename T> inline Vector DynLanczos<Matrix, T>::getUpdatedEigenvector(int i) {
        this->assureFinished();
        assert(i < basis.size());
        return basis[i];
        
  }
    

    
    template<class Matrix, typename T> inline void DynLanczos<Matrix, T>:: loadEigenvectors(std::vector<Vector> in_eigenvectors) {
        this->assureFinished();
        const int n = Lanczos<Matrix,T>::A.numberOfRows();
        assert(basisOld.size() == in_eigenvectors.size());
        for (int i = 0; i < in_eigenvectors.size(); i++) {
            if (basisOld[i].getDimension() != n)
                basisOld[i] = Vector(in_eigenvectors[i]);
            else
                basisOld[i] = in_eigenvectors[i];
        }
    }
        
  

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_DYN_LANCZOS_HPP_
