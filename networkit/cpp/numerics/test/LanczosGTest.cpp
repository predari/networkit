#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/Eigen.hpp>

namespace NetworKit {

class LanczosGTest : public testing::Test {};

TEST(LanczosGTest, debugkEigen) {
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
    
    Lanczos<CSRMatrix,double> s(A, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    // get eigenvectors without computing them first.
    evectors = s.getkEigenvectors();
    INFO("CHECK!!  evectors_size (not computed yet):" , evectors.size());
    for (int i =0; i< evectors.size(); i++) {
      INFO("CHECK evectors[,", i, ",]_size :", evectors[i].getDimension());
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    INFO("CHECK!!  evectors_size (computed!!):" , evectors.size());
    for (int i =0; i< evectors.size(); i++) {
      INFO("vec ", i);
      for (int j = 0; j < evectors[i].getDimension(); j++) {  
	INFO(evectors[i][j]);
      }
    } 

    
	
    Vector vec(Vector(4,0.));
    vec = s.getEigenvector(0);

    INFO("Eigenvalues of matrix A in R4x4:");
    for (int i =0; i< k; i++)  
      INFO("e[",i,"] = ", e[i]);
    
    // INFO("Eigenvector of matrix A evector[0] :", vec.getDimension());
    // ASSERT_EQ(4u, vec.getDimension());
    // for (int i =0; i< vec.getDimension(); i++) {  
    //   INFO(" (vec vs v) : ", vec[i] , ",", v[0][i] );
    //   //EXPECT_NEAR(vec[i], v[0][i], 1e-5);
    // }

    ASSERT_EQ(k, evectors.size());
    INFO(" b size = ", evectors.size());
    for (int i = 0; i < evectors.size(); i++) {
      ASSERT_EQ(n, evectors[i].getDimension());
      INFO(" k = ", i);
      ASSERT_EQ(n, evectors[i].getDimension());
      for (int j = 0; j < evectors[i].getDimension(); j++) {  
	INFO(" (vec vs real) : ", evectors[i][j] , "  ,  ", v[i][j] );
    	//EXPECT_NEAR(evectors[i][j], v[i][j], 1e-5);
      }
    }
    
    for (int i = 0; i < eigens.size(); i++) {
      EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    if (!s.checkEigenvectors())
      std::cout << "*** WHY?????? ***" << std::endl;
    

}


TEST(LanczosGTest, debugkEigenGraph) {
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

    const auto L = CSRMatrix::laplacianMatrix(G);
    
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    // actual eigenvalues and eigenvectors

         
    std::vector<Vector> evectors(k,Vector(4,0.));
    
    Lanczos<CSRMatrix,double> s(L, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    // get eigenvectors without computing them first.
    evectors = s.getkEigenvectors();
    INFO("CHECK!!  evectors_size (not computed yet):" , evectors.size());
    for (int i =0; i< evectors.size(); i++) {
      INFO("CHECK evectors[,", i, ",]_size :", evectors[i].getDimension());
    } 
    s.computekEigenvectors();
    evectors = s.getkEigenvectors();
    INFO("CHECK!!  evectors_size (computed!!):" , evectors.size());
    for (int i =0; i< evectors.size(); i++) {
      INFO("vec ", i);
      for (int j = 0; j < evectors[i].getDimension(); j++) {  
	INFO(evectors[i][j]);
      }
    } 

    
	
    Vector vec(Vector(4,0.));
    vec = s.getEigenvector(0);

    INFO("Eigenvalues of matrix A in R4x4:");
    for (int i =0; i< k; i++)  
      INFO("e[",i,"] = ", e[i]);
    

    // ASSERT_EQ(k, evectors.size());
    // INFO(" b size = ", evectors.size());
    // for (int i = 0; i < evectors.size(); i++) {
    //   ASSERT_EQ(n, evectors[i].getDimension());
    //   INFO(" k = ", i);
    //   ASSERT_EQ(n, evectors[i].getDimension());
    //   for (int j = 0; j < evectors[i].getDimension(); j++) {  
    // 	INFO(" (vec vs real) : ", evectors[i][j] , "  ,  ", v[i][j] );
    // 	//EXPECT_NEAR(evectors[i][j], v[i][j], 1e-5);
    //   }
    // }
    
    // for (int i = 0; i < eigens.size(); i++) {
    //   EXPECT_NEAR(eigens[i], e[i], 1e-5);
    // }

}




  
TEST(LanczosGTest, debugkEigenDynamic) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8

    
    DynamicMatrix A(4, triplets);
    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);

    std::vector<double> eigens = {14.0734777528174, 10.8190609243708, 8.14343524715178, 5.96402607566004 };
    const std::vector<std::vector<double>> v = {{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
					        { 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
						{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 }};

    std::vector<Vector> b(k,Vector(4,0.));
    
    
    Lanczos<DynamicMatrix,double> s(A, k, a, skip, true);
    s.run();
    e = s.getkEigenvalues();
    // get eigenvectors without computing them first.
    b = s.getkEigenvectors();
    

    
    EXPECT_NEAR(eigens[0], e[0], 1e-5);
    EXPECT_NEAR(eigens[1], e[1], 1e-5);

}

//   // COMMENT: Does not work because of Vector implemented only for double. 
// TEST(LanczosGTest, debugkEigenSinglePrecision) {
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
//     // TODO: compute proper values for result
//     std::vector<float> result = {6, 13};
//     std::vector<Vector> b(k,Vector(4,0.));
    
//     Lanczos<CSRMatrix,float> s(A, k, a, skip);
//     //s.run();
//     //e = s.getAllEigenvalues();
//     //b = s.getBasis();
    
//     //EXPECT_NEAR(6.0, e[0], 1e-5);
//     //EXPECT_NEAR(13.0, e[1], 1e-5);
// }


  
TEST(LanczosGTest, debugFullEigen) {
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
    const std::vector<std::vector<double>> v = {{ 0.393474513266299, -0.680207701966034, 0.461300987062449, -0.411942579652973 },
						{-0.623361112531999, -0.491294682496488, -0.500932705771864, -0.345133137530461},
					        { 0.640429974959434, -0.226704300124538, -0.707852009150845, 0.193391159621156 },
						{-0.215518393568104, -0.494522170829382, 0.187652222913562, 0.820844862216739 }};
    
    std::vector<Vector> b(n,Vector(n,0.));
    Lanczos<CSRMatrix,double> s(A, n, a, skip, true);
    s.run();
    // e = s.getAllEigenvalues();
    // b = s.getBasis();


    // INFO("Eigenvalues of matrix A in R4x4:");
    // for (int i =0; i< n; i++)  
    //   INFO("e[",i,"] = ", e[i]);
    
    // INFO("Eigenvectors of matrix A in R4x4:");
    // for (int i =0; i< n; i++) {  
    //   INFO("vec[",i,"] : ");
    //   Vector r = s.getEigenvector(i);
    //   if (r != b[i])
    // 	ERROR(" Eigen::getEigenvector() does not work as intended. ");
    //   int siz = b[i].length();
    //   if (siz != A.numberOfColumns())
    // 	ERROR(" Problem with size of eigenvector! Not equal to number of cols of A.");
    //   for (int j =0; j < siz; j ++) {
    // 	//INFO(r[j]);
    // 	EXPECT_NEAR(r[j], v[i][j], 1e-5);

    //   }
    // }

    
    EXPECT_NEAR(eigens[0], e[0], 1e-5);
    EXPECT_NEAR(eigens[1], e[1], 1e-5);
    EXPECT_NEAR(eigens[2], e[2], 1e-5);
    EXPECT_NEAR(eigens[3], e[3], 1e-5);

}






  
} /* namespace NetworKit */
