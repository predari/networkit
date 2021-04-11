#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/Lanczos.hpp>

namespace NetworKit {

class LanczosGTest : public testing::Test {};



//     std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
// //	10  -1   2   0
// //	-1  11  -1   3
// //	 2  -1  10  -1
// //	 0   3  -1   8
//     // actual eigenvalues and eigenvectors
//     std::vector<double> eigens = {14.0734777528174, 10.8190609243708 };
//     const std::vector<std::vector<double>> v = {{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 },
// 						{ 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
// 						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
// 						{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
// 						};

    
TEST(LanczosGTest, kEigenvectorsMatrix) {
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
    ASSERT_EQ(k, e.size());
    ASSERT_LE(e.size(), n);
    
    // for (int i = 0; i < e.size(); i++) {
    //   EXPECT_NEAR(eigens[i], e[i], 1e-5);
    // }
    std::vector<Vector> evectors(k,Vector(n,0.));
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    if (!s.checkEigenvectors())
      WARN(" Eigenvectors are not correct!");
    
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

    

    // actual eigenvalues and eigenvectors
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

    int j = 1;
    ASSERT_LE(j, k);
    Vector vec(Vector(n,0.));
    vec = s.getEigenvector(j);
    ASSERT_EQ(n, vec.getDimension());
    for (int i = 0; i< vec.getDimension(); i++) {  
      INFO(" (comp vs v) : ", vec[i] , ",", v2[i] );
      EXPECT_NEAR(fabs(vec[i]), fabs(v2[i]), 1e-5);
    }
    
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

    

    // actual eigenvalues and eigenvectors
    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };
    const std::vector<std::vector<double>> v = { {-0.195439507584855,-0.195439507584855, 0.827895039618530, -0.316227766016838, -0.316227766016838, 0.195439507584855},
                                                 {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221},
                                                 {-4.82186041151648e-16, -4.82186041151648e-16,3.33066907387547e-16, -0.707106781186547, 0.707106781186547, -8.32667268468868e-16},
                                                 { 0.707106781186548, -0.707106781186548, 6.16297582203916e-32, 7.85046229341889e-17, -7.85046229341889e-17, 0},
                                                 { 0.511667273601693, 0.511667273601693, 0.120788258431983, -0.316227766016838, -0.316227766016838, -0.511667273601693},
                                                 { 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863 }};


    

    const auto L = CSRMatrix::laplacianMatrix(G);

    INFO("L[2,2] = ", L(2,2));
    INFO("L[0,2] = ", L(0,2));
    INFO("L[2,0] = ", L(2,0));
    INFO("L[4,0] = ", L(4,0));
    INFO("L[0,4] = ", L(0,4));
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    
    Lanczos<CSRMatrix,double> s(L, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    ASSERT_EQ(k, e.size());
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      INFO(e[i]);
    }
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    evectors = s.getkEigenvectors();
    for (int i =0; i< evectors.size(); i++) {
        if (!evectors[i].getDimension()) {
            WARN(" Eigenvectors are not yet computed. Please call computeEigenvectors() first.");
            break;
        }
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    if (!s.checkEigenvectors())
      WARN(" Eigenvectors are not correct! :(");
    ASSERT_TRUE(s.checkEigenvectors());
    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        INFO(" eigenvector : ", i );
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            INFO(" (comp vs v) : ", evectors[i][j] , ",", v[i][j] );
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-6);
      }
    }    
}




TEST(LanczosGTest,  kEigenvaluesDynamicLaplacian) {
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

    

    // actual eigenvalues and eigenvectors
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
    ASSERT_EQ(k, e.size());
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      INFO(e[i]);
    }
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    evectors = s.getkEigenvectors();
    for (int i =0; i< evectors.size(); i++) {
        if (!evectors[i].getDimension()) {
            WARN(" Eigenvectors are not yet computed. Please call computeEigenvectors() first.");
            break;
        }
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    if (!s.checkEigenvectors())
      WARN(" Eigenvectors are not correct! :(");
    ASSERT_TRUE(s.checkEigenvectors());
    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        INFO(" eigenvector : ", i );
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            INFO(" (comp vs v) : ", evectors[i][j] , ",", v[i][j] );
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-6);
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
    ASSERT_EQ(n, e.size());
    ASSERT_LE(e.size(), eigens.size());
    
    for (int i = 0; i < e.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

}






  
} /* namespace NetworKit */
