#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/DynLanczosEigen.hpp>
#include <networkit/auxiliary/Random.hpp>
//#include <networkit/algebraic/AlgebraicGlobals.hpp>

//#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class DynLanczosGTest : public testing::Test {};



TEST(DynLanczosGTest, runDynVsStatic) {

  std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {0,6,1}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {1,5,6}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {2,4,3}, {2,5,7}, {3,1,3}, {3,2,-1}, {3,3,8}, {3,4,2}, {4,0,2}, {5,1,4}, {5,5,5}, {6,3,6}, {6,6,9}};
//	10  -1   2   0  0  0  1
//	-1  11  -1   3  0  6  0
//	 2  -1  10  -1  3  7  0
//	 0   3  -1   8  2  0  0
//       1   0   0   0  0  0  0
//       0   4   0   0  0  5  0
//       0   0   0   6  0  0  9  
  int n = 7;
  CSRMatrix A(n, triplets);

  int k = 3;
  int a = 2;
  int skip = 16;
  std::vector<double> e(k,0.);

  DEBUG("Initializing DynLanczos");
  DynLanczos<CSRMatrix,double> dynl(A, k, a, skip);
  DEBUG("Initializing Lanczos");
  Lanczos<CSRMatrix, double> l(A, k, a, skip);

  DEBUG("Running DynLanczos");
  dynl.run();
  DEBUG("Running Lanczos");
  l.run();
  std::vector<double> dynl_eigen = dynl.getkEigenvalues();
  std::vector<double> l_eigen = l.getkEigenvalues();
    double err1=0;
    for(count i=0; i<k; i++) {
        double x = dynl_eigen[i]-l_eigen[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG(" err1 = ", err1);
    DEBUG("Before the edge insertion: ");
    std::vector<GraphEvent> batch;
    count nInsertions = 5, i = 0;
    while (i < nInsertions) {
      uint64_t v1 = Aux::Random::integer(0,n-1);
      uint64_t v2 = Aux::Random::integer(0,n-1);
      if (v1 != v2 && !A(v1, v2)) {
	double weight = 1.0;
	DEBUG("Adding element ", v1, ", ", v2, " Before: m[v1,v2] = ", A(v1, v2));
	A.setValue(v1, v2, weight);
	//A.setValue(v2, v1, weight);
	DEBUG(" After: m[v1,v2] = ", A(v1, v2));
	batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, weight));
	DEBUG(" After push.back");
	i++;
      }
    }
    DEBUG("Running Lanczos (again)");
    l.run();
    DEBUG("Updating DynLanczos");
    dynl.updateBatch(batch);
    DEBUG("Calling DynLanczos Eigenvalues");
    dynl_eigen = dynl.getkEigenvalues();
    DEBUG("Calling nLanczos Eigenvalues");
    l_eigen = dynl.getkEigenvalues();
    err1 = 0;
    for(count i=0; i<k; i++) {
      double x = dynl_eigen[i]-l_eigen[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG("After the edge insertion: ");
}


  


  
} /* namespace NetworKit */
