#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/Eigen.hpp>

namespace NetworKit {

class LanczosGTest : public testing::Test {};

TEST(LanczosGTest, kEigenvaluesMatrix) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8
    CSRMatrix A(4, triplets);
    int n = A.numberOfRows();
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    // actual eigenvalues and eigenvectors
    std::vector<double> eigens = {14.0734777528174, 10.8190609243708 };
    Lanczos<CSRMatrix,double> s(A, k, a, skip);
    s.run();
    e = s.getkEigenvalues();
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

}

TEST(LanczosGTest, OneEigenvectorMatrix) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8
    CSRMatrix A(4, triplets);
    int n = A.numberOfRows();
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    // actual eigenvalues and eigenvectors
    std::vector<double> eigens = {14.0734777528174, 10.8190609243708 };
    const std::vector<std::vector<double>> v = {{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 },
						{ 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
						{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
						};
         
    Lanczos<CSRMatrix,double> s(A, k, a, skip, true, 1e-8);
    s.run();
    e = s.getkEigenvalues();
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    int j = 1;
    Vector vec(Vector(4,0.));
    vec = s.getEigenvector(eigens[j]);

    INFO("Eigenvector of matrix A evector[0] :", vec.getDimension());
    ASSERT_EQ(4u, vec.getDimension());    
    for (int i = 0; i< vec.getDimension(); i++) {  
      INFO(" (comp vs v) : ", vec[i] , ",", v[j][i] );
      //EXPECT_NEAR(vec[i], v[0][i], 1e-8);
    }
    
}
  

TEST(LanczosGTest, kEigenvectorsMatrix) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8
    CSRMatrix A(4, triplets);
    int n = A.numberOfRows();
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    // actual eigenvalues and eigenvectors
    std::vector<double> eigens = {14.0734777528174, 10.8190609243708 };
    const std::vector<std::vector<double>> v = {{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 },
						{ 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
						{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
						};

         
    std::vector<Vector> evectors(k,Vector(4,0.));
    Lanczos<CSRMatrix,double> s(A, k, a, skip);
    s.run();
    e = s.getkEigenvalues();
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }


    evectors = s.getkEigenvectors();
    for (int i =0; i< evectors.size(); i++) {
        if (!evectors[i].getDimension()) {
            WARN(" Eigenvectors are not yet computed. Please call computeEigenvectos()");
            break;
        }
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    if (!s.checkEigenvectors())
      WARN(" Eigenvectors are not correct!");
    
    for (int i =0; i< evectors.size(); i++) {
      ASSERT_EQ(n, evectors[i].getDimension());
      for (int j = 0; j< evectors[i].getDimension(); j++) {
          INFO(" (comp vs v) : ", evectors[i][j] , ",", v[i][j] );
          //EXPECT_NEAR(evectors[i][j], v[i][i], 1e-8);
      }
    }
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


    // G.addEdge(0, 2);
    // G.addEdge(1, 2);
    // G.addEdge(2, 3);
    // G.addEdge(2, 4);
    // G.addEdge(3, 5);
    // G.addEdge(4, 5);

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
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    int j = 2;
    Vector vec(Vector(n,0.));
    vec = s.getEigenvector(eigens[j]);
    ASSERT_EQ(n, vec.getDimension());
    for (int i = 0; i< vec.getDimension(); i++) {  
      INFO(" (comp vs v) : ", vec[i] , ",", v2[i] );
      //EXPECT_NEAR(vec[i], v[0][i], 1e-8);
    }
    
}




  
TEST(LanczosGTest, kEigenvaluesDynamicMatrix) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8

    DynamicMatrix A(4, triplets);
    int n = A.numberOfRows();
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    // actual eigenvalues and eigenvectors
    std::vector<double> eigens = {14.0734777528174, 10.8190609243708 };
    const std::vector<std::vector<double>> v = {{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 },
						{ 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
						{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
						};

         
    std::vector<Vector> evectors(k,Vector(4,0.));
    Lanczos<DynamicMatrix,double> s(A, k, a, skip);
    s.run();
    e = s.getkEigenvalues();
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }


    evectors = s.getkEigenvectors();
    for (int i =0; i< evectors.size(); i++) {
      if (!evectors[i].getDimension()) 
	WARN(" Eigenvectors are not yet computed. Please call computeEigenvectos()");
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    if (!s.checkEigenvectors())
      WARN(" Eigenvectors are not correct!");

    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            INFO(" (comp vs v) : ", evectors[i][j] , ",", v[i][j] );
            //EXPECT_NEAR(evectors[i][j], v[i][i], 1e-8);
      }
    }
}


//   // COMMENT: Does not work because of Vector implemented only for double. 
// TEST(LanczosGTest, kEigenValuesSinglePrecision) {
//     std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
// //	10  -1   2   0
// //	-1  11  -1   3
// //	 2  -1  10  -1
// //	 0   3  -1   8
//     CSRMatrix A(4, triplets);
//     int k = 2;
//     int a = 2;
//     int skip = 16;
//     std::vector<float> e(k,0.);
//     std::vector<float> eigens = {14.0734777528174, 10.8190609243708, 8.14343524715178, 5.96402607566004 };
    
//     Lanczos<CSRMatrix,float> s(A, 3, a, skip, false);
//     s.run();
//     e = s.getkEigenvalues();
//     for (int i = 0; i < eigens.size(); i++) {
//       EXPECT_NEAR(eigens[i], e[i], 1e-5);
//     }
    
// }


  
TEST(LanczosGTest, FullEigenvaluesMatrix) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8
    CSRMatrix A(4, triplets);
    int n = A.numberOfRows();
    int a = 2;
    int skip = 16;
    std::vector<double> e(n,0.);

    std::vector<double> eigens = {14.0734777528174, 10.8190609243708, 8.14343524715178, 5.96402607566004 };
    
    std::vector<Vector> b(n,Vector(n,0.));
    Lanczos<CSRMatrix,double> s(A, n, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }
}






  
} /* namespace NetworKit */
