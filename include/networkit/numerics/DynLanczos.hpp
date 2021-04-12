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
 * @ingroup centrality
 * Dynamic APSP.
 */

template<class Matrix, typename T>
class DynLanczos : public Lanczos<Matrix, T>, public DynAlgorithm {

  
public:
    /**
     * Creates the object for @a G.
     *
     * @param G The graph.
     */
    //DynLanczos(const Matrix & A) :
    //  Lanczos<Matrix,double>(A) {}
    DynLanczos(const Matrix & A, const int k, int a, int skip, bool double_precision = true);

    ~DynLanczos() override = default;

    /**
     * Runs static betweenness centrality algorithm on the initial graph.
     */
    void run() override;
    
    /**
     * Updates the betweenness centralities after an edge insertions on the graph.
     * Notice: it works only with edge insertions.
     *
     * @param e The edge insertions.
     */
    void update(GraphEvent singleEvent) override {
        std::vector<GraphEvent> events{singleEvent};
        updateBatch(events);
    }
    /**
     * Updates the betweenness centralities after a batch of edge insertions on the graph.
     * Notice: it works only with edge insertions.
     *
     * @param batch The batch of edge insertions.
     */
    
    void updateBatch(const std::vector<GraphEvent>& batch) override;

   

    std::vector<T> getAllEigenvalues() ;
    
    T getEigenvalue(int i) ;
    
    std::vector<Vector> getBasis() ;
  
    Vector getEigenvector(int i) ;

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
  DynLanczos<Matrix, T>::DynLanczos(const Matrix & A, const int k, int a, int skip, bool double_precision) :
      Lanczos<Matrix,T>(A, k, a, skip, double_precision), Delta(Matrix(A.numberOfRows(), A.numberOfColumns())), eigen(k, 0.), basis(k, Vector(A.numberOfRows(), 0.)){ }
   

  
template<class Matrix, typename T> void DynLanczos<Matrix, T> :: run() {

    Lanczos<Matrix,T>::run();
    Lanczos<Matrix,T>::assureFinished();
    Lanczos<Matrix,T>::computekEigenvectors();
    
   
    eigenOld = Lanczos<Matrix,T>::eigenvalues;
    basisOld = Lanczos<Matrix,T>::eigenvectors;
    int n = Lanczos<Matrix,T>::A.numberOfColumns();
    int k = Lanczos<Matrix,T>::k;
    assert(eigenOld.size() == k);
    assert(basisOld.size() == k);

    DEBUG("[DYN] ************************************** ");
    DEBUG("[DYN] ************* run() ****************** ");
    std::cout << "[ ";
    for (int i = 0; i < eigenOld.size(); i++)
        std::cout << eigenOld[i] << " ";
    std::cout << "]\n";
    DEBUG("[DYN] ************************************** ");
    for (int i =0; i< basisOld.size(); i++) {
        assert(basisOld[i].getDimension() == n);
        std::cout << " [DYN] ******** Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basisOld[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("[DYN] ************************************** ");


    
    //eigenOld = Lanczos::getEigenvalues();
    //basisOld = Lanczos::getBasis();
    this->hasRun = true;
}


  

  // TODOs:
  // 1. change double for eigenvalues with template T
  // 2. update the previous results to be the old ones
  // 3. allow edge insertions/deletions and edge weight decreases and increases
  // 4. change private variables to protected
  // 5. consider parallelism
  template<class Matrix, typename T> void DynLanczos<Matrix, T> :: updateBatch(const std::vector<GraphEvent>& batch) {
    INFO("Entering update");
    // create a delta-update matrix
    const int n = Lanczos<Matrix,T>::A.numberOfRows();
    INFO(" size of original matrix n       = ", n);
    INFO(" nnz of original matrix          = ", Lanczos<Matrix,T>::A.nnz());
    const int batchSize = batch.size();

    INFO(" size of delta matrix n (before) = ", Delta.numberOfRows());
    INFO(" nnz of delta matrix (before)    = ", Delta.nnz());
    assert(!Delta.nnz());
   
    
    for(auto event : batch){
      node u = event.u;
      node v = event.v;
      if (!(event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
      }
            
      Delta.setValue(u, v, -1.0);
      Delta.setValue(v, u, -1.0);
      Delta.setValue(v, v, Delta(v, v) + event.w);
      Delta.setValue(u, u, Delta(u, u) + event.w);
    }

    INFO(" size of delta matrix n (after)  = ", Delta.numberOfRows());
    INFO(" nnz of delta matrix (after)     = ", Delta.nnz());
    
    //for(int i = 0; i < triplets.size(); i++)
    //INFO(" edge ", triplets[i].row , ",", triplets[i].column);
    //INFO(" Delta.rows() = ", Delta.numberOfRows() , " Delta.cols() = ", Delta.numberOfColumns());
    //INFO(" Delta.nnz() = ", Delta.nnz());


    // INFO("D[1,3] = ", Delta(1,3));
    // INFO("D[3,3] = ", Delta(3,3));
    // INFO("D[4,4] = ", Delta(4,4));
    // INFO("D[1,1] = ", Delta(1,1));
    // INFO("D[1,2] = ", Delta(1,2));

    
    for(int j = 0; j < Lanczos<Matrix,T>::k; j++) {
      Vector deltaEigenvector(n,0.);
      Vector v(n,0.);
      for(int i = 0; i < Lanczos<Matrix,T>::k; i++) {
          if (i == j) continue;
      
          Vector b = basisOld[j];
          assert(b.getDimension() == n);
          
          v = Delta * basisOld[j];
          double alpha = basisOld[i].transpose() * v;
          alpha /=  (eigenOld[j] - eigenOld[i]);
          INFO(" delta eigenvector");
          deltaEigenvector = deltaEigenvector + alpha * basisOld[i];
          for (int i= 0; i < n; i++)
              INFO(deltaEigenvector[i]);
          
      }
      double deltaEigenvalue = basisOld[j].transpose() * v;
      eigen[j] = eigenOld[j] + deltaEigenvalue; 
      basis[j] = basisOld[j] + deltaEigenvector;
    }


    DEBUG("[DYN] [OLD] ******************************** ");
    DEBUG("[DYN] [OLD] ********* update() ************* ");
    std::cout << "[ ";
    for (int i = 0; i < eigenOld.size(); i++)
        std::cout << eigenOld[i] << " ";
    std::cout << "]\n";
    DEBUG("[DYN] [OLD] ******************************** ");
    for (int i =0; i< basisOld.size(); i++) {
        assert(basisOld[i].getDimension() == n);
        std::cout << " [DYN] [OLD] *** Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basisOld[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("[DYN] ************************************** ");

    // TODO: change as to receive Delta as input argument in update
    // update for next run
    // clear Delta
    eigenOld = eigen;
    basisOld = basis;
    Matrix M = Matrix(Delta.numberOfRows());
    Lanczos<Matrix,T>::setup(Lanczos<Matrix,T>::A + Delta);
    Delta = M;


    DEBUG("[DYN] [NEW] ******************************** ");
    DEBUG("[DYN] [NEW] ********* update() ************* ");
    std::cout << "[ ";
    for (int i = 0; i < eigen.size(); i++)
        std::cout << eigen[i] << " ";
    std::cout << "]\n";
    DEBUG("[DYN] [NEW] ******************************** ");
    for (int i =0; i< basis.size(); i++) {
        assert(basis[i].getDimension() == n);
        std::cout << " [DYN] [NEW]  *** Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << basis[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("[DYN] ************************************** ");


    
  }


  template<class Matrix, typename T> inline std::vector<T> DynLanczos<Matrix, T> :: getAllEigenvalues() {
    this->assureFinished();
    return eigen;
}

  template<class Matrix, typename T> inline T DynLanczos<Matrix, T>::getEigenvalue(int i) {
    this->assureFinished();
    return eigen[i];
}

  template<class Matrix, typename T> inline  std::vector<Vector> DynLanczos<Matrix, T>::getBasis() {
    this->assureFinished();
    return basis;
}

    
  template<class Matrix, typename T> inline Vector DynLanczos<Matrix, T>::getEigenvector(int i) {
    this->assureFinished();
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
