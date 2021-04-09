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
#include <networkit/numerics/Eigen.hpp>
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
  DynLanczos(const Matrix & A, const int k, int a, int skip, bool double_precision = false);

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
  // COMMENT: why not override to override those of Lanczos ? 
  void updateBatch(const std::vector<GraphEvent>& batch);

  std::vector<T> getAllEigenvalues() ;
    
  T getEigenvalue(int i) ;
  
  std::vector<Vector> getBasis() ;
  
  Vector getEigenvector(int i) ;
  

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
    Lanczos<Matrix,T>(A, k, a, skip, double_precision), Delta(Matrix()), eigen(k, 0.), basis(k, Vector(A.numberOfRows(), 0.)){ }
   

  
template<class Matrix, typename T> void DynLanczos<Matrix, T> :: run() {
  Lanczos<Matrix,T>::run();
  eigenOld = Lanczos<Matrix,T>::eigen;
  basisOld = Lanczos<Matrix,T>::basis; 
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
    INFO(" size of original matrix n = ", n);
    const int batchSize = batch.size();
    std::vector<Triplet> triplets;

    
    for(auto event : batch){
      node u = event.u;
      node v = event.v;
      //double entry = A(u,v);
      if (!(event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
      }
      triplets.push_back({u, v, event.w});
    }
    for(int i = 0; i < triplets.size(); i++)
      INFO(" edge ", triplets[i].row , ",", triplets[i].column);

    INFO(" Delta.rows() = ", Delta.numberOfRows() , " Delta.cols() = ", Delta.numberOfColumns());
    INFO(" Delta.nnz() = ", Delta.nnz());

    Matrix m(n,triplets);
    Delta = m;
    
    INFO(" Delta.rows() = ", Delta.numberOfRows() , " Delta.cols() = ", Delta.numberOfColumns());
    INFO(" Delta.nnz() = ", Delta.nnz());
    for(int i = 0; i < Delta.numberOfRows(); i++)
      INFO(" Delta.nnzPerRow() = ", Delta.nnzInRow(i));

    
    INFO(" k = ", Lanczos<Matrix,T>::k);
    
    for(int j = 0; j < Lanczos<Matrix,T>::k; j++) {
      Vector deltaEigenvector(n,0.);
      Vector v(n,0.);
      for(int i = 0; i < Lanczos<Matrix,T>::k; i++) {
	if (i == j) continue;
	INFO(" -- Basis[j]:");

	Vector b = basisOld[j];
	INFO(" -- lenght = ", b.length());	
	for (int i= 0; i < n; i++)
	  INFO(b[i]);
	
	v = Delta * basisOld[j];
	INFO(" -- A*v");
	double alpha = basisOld[i].transpose() * v;
	INFO(" -- v^T*v");
	alpha /= (eigenOld[j] - eigenOld[i]);
	INFO(" alpha");
        deltaEigenvector = deltaEigenvector + alpha * basisOld[i]; 
      }
      double deltaEigenvalue = basisOld[j].transpose() * v;
      eigen[j] = eigenOld[j] + deltaEigenvalue; 
      basis[j] = basisOld[j] + deltaEigenvector;
      	INFO(" delta updates");
    }

    eigenOld = eigen;
    basisOld = basis;    
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


  

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_DYN_LANCZOS_HPP_
