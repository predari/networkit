#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/Lanczos.hpp>

namespace NetworKit {

class LanczosGTest : public testing::Test {};


    
TEST(LanczosGTest, kEigenvectorsAdjacency) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    const auto A = CSRMatrix::adjacencyMatrix(G);    
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);


    Lanczos<CSRMatrix,double> s(A, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_LE(e.size(), k);
    
    s.computekEigenvectors();
    std::vector<Vector> evectors = s.getkEigenvectors();
    ASSERT_EQ(evectors.size(),e.size());
    ASSERT_TRUE(s.checkEigenvectors());
    
    
}


  

TEST(LanczosGTest,  OneEigenvectorLaplacian) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();


    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    
    // eigenpairs computed with matlab
    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };
    const std::vector<double> v2 = {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221 };
    

    const auto L = CSRMatrix::laplacianMatrix(G);

    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
         
    std::vector<Vector> evectors(k,Vector(n,0.));
    
    Lanczos<CSRMatrix,double> s(L, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      INFO(e[i]);
    }
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    Vector vec(Vector(n,0.));
    // get first eigenvector
    vec = s.getEigenvector(1);
    ASSERT_EQ(n, vec.getDimension());
    for (int i = 0; i< vec.getDimension(); i++)
      EXPECT_NEAR(fabs(vec[i]), fabs(v2[i]), 1e-5);
    
}


TEST(LanczosGTest,  kEigenvectorsLaplacian) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();


    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    

    // eigenpairs computed with matlab
    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };
    const std::vector<std::vector<double>> v = { {-0.195439507584855,-0.195439507584855, 0.827895039618530, -0.316227766016838, -0.316227766016838, 0.195439507584855},
                                                 {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221},
                                                 {-4.82186041151648e-16, -4.82186041151648e-16,3.33066907387547e-16, -0.707106781186547, 0.707106781186547, -8.32667268468868e-16},
                                                 { 0.707106781186548, -0.707106781186548, 6.16297582203916e-32, 7.85046229341889e-17, -7.85046229341889e-17, 0},
                                                 { 0.511667273601693, 0.511667273601693, 0.120788258431983, -0.316227766016838, -0.316227766016838, -0.511667273601693},
                                                 { 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863 }};


    

    const auto L = CSRMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    
    Lanczos<CSRMatrix,double> s(L, k, a, skip, true, 1e-06);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      INFO(e[i]);
    }
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    ASSERT_EQ(evectors.size(),e.size());
    ASSERT_TRUE(s.checkEigenvectors());
    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-5);
      }
    }    
}




TEST(LanczosGTest,  kEigenvectorsDynamicLaplacian) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();


    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    

    // eigenpairs computed with matlab
    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };
    const std::vector<std::vector<double>> v = { {-0.195439507584855,-0.195439507584855, 0.827895039618530, -0.316227766016838, -0.316227766016838, 0.195439507584855},
                                                 {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221},
                                                 {-4.82186041151648e-16, -4.82186041151648e-16,3.33066907387547e-16, -0.707106781186547, 0.707106781186547, -8.32667268468868e-16},
                                                 { 0.707106781186548, -0.707106781186548, 6.16297582203916e-32, 7.85046229341889e-17, -7.85046229341889e-17, 0},
                                                 { 0.511667273601693, 0.511667273601693, 0.120788258431983, -0.316227766016838, -0.316227766016838, -0.511667273601693},
                                                 { 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863 }};


    

    const DynamicMatrix L = DynamicMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    
    Lanczos<DynamicMatrix,double> s(L, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      INFO(e[i]);
    }
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    ASSERT_EQ(evectors.size(),e.size());
    ASSERT_TRUE(s.checkEigenvectors());
    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-5);
      }
    }

    
}


  
TEST(LanczosGTest, FullEigenvaluesLaplacian) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();


    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    

    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };

    

    const auto L = CSRMatrix::laplacianMatrix(G);
    int a = 2;
    int skip = 16;
    std::vector<double> e(n,0.);
    
    Lanczos<CSRMatrix,double> s(L, n, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

}






  
} /* namespace NetworKit */
