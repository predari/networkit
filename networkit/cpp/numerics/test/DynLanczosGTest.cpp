#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/DynLanczos.hpp>
#include <networkit/auxiliary/Random.hpp>
//#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/graph/GraphTools.hpp>

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


TEST(DynLanczosGTest, DynamicVsStaticOnStaticGraph) {

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

    // eigenvalues and eigenvectors computed with matlab
    std::vector<double> eigens = {5.2360679, 3.000, 2.0000, 1.0000, 0.763932022500210, 2.83534217966776e-17 };
    const std::vector<std::vector<double>> v = { {-0.195439507584855,-0.195439507584855, 0.827895039618530, -0.316227766016838, -0.316227766016838, 0.195439507584855},
                                                 {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221},
                                                 {-4.82186041151648e-16, -4.82186041151648e-16,3.33066907387547e-16, -0.707106781186547, 0.707106781186547, -8.32667268468868e-16},
                                                 { 0.707106781186548, -0.707106781186548, 6.16297582203916e-32, 7.85046229341889e-17, -7.85046229341889e-17, 0},
                                                 { 0.511667273601693, 0.511667273601693, 0.120788258431983, -0.316227766016838, -0.316227766016838, -0.511667273601693},
                                                 { 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863 }};


    

    auto L = CSRMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    
    // run static and compare with matlab 
    Lanczos<CSRMatrix, double> l(L, k, a, skip, true, 1e-05);
    l.run();
    e = l.getkEigenvalues();
    ASSERT_LE(e.size(), k);
    for (int i = 0; i < e.size(); i++)
        INFO(e[i]);
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    ASSERT_TRUE(l.checkEigenvectors());
    for (int i = 0; i < evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-5);
        }
    }    
    
    // run dynamic and compare with matlab (no update)
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip, true, 1e-05);
    dyn_l.run();
    std::vector<double> dyn_e = dyn_l.getkEigenvalues();
    ASSERT_LE(dyn_e.size(), k);
    
    double err1= 0;
    for(count i=0; i<k; i++) {
        double x = dyn_e[i]-e[i];
        if (x > err1)
            err1 = x;
    }
    if (err1)
        WARN(" Static and dynamic should produce the same result on a static graph.");
    
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(eigens[i], dyn_e[i], 1e-5);
        EXPECT_NEAR(e[i], dyn_e[i], 1e-5);
    }

    std::vector<Vector> dyn_evectors(k,Vector(n,0.));
    dyn_l.computekEigenvectors();
    dyn_evectors = dyn_l.getkEigenvectors();
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(v[i][j]), 1e-5);
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(evectors[i][j]), 1e-5);
        }
    }    
}

    


TEST(DynLanczosGTest, runDynVsStaticLaplacian) {

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
    const std::vector<std::vector<double>> v = { {-0.195439507584855,-0.195439507584855, 0.827895039618530, -0.316227766016838, -0.316227766016838, 0.195439507584855},
                                                 {0.182574185835055, 0.182574185835055, -0.365148371670111, -0.365148371670111, -0.365148371670110, 0.730296743340221},
                                                 {-4.82186041151648e-16, -4.82186041151648e-16,3.33066907387547e-16, -0.707106781186547, 0.707106781186547, -8.32667268468868e-16},
                                                 { 0.707106781186548, -0.707106781186548, 6.16297582203916e-32, 7.85046229341889e-17, -7.85046229341889e-17, 0},
                                                 { 0.511667273601693, 0.511667273601693, 0.120788258431983, -0.316227766016838, -0.316227766016838, -0.511667273601693},
                                                 { 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863, 0.408248290463863 }};


    

    auto L = CSRMatrix::laplacianMatrix(G);

    //INFO("L[2,2] = ", L(2,2));
    //INFO("L[0,2] = ", L(0,2));
    //INFO("L[2,0] = ", L(2,0));
    //INFO("L[4,0] = ", L(4,0));
    //INFO("L[0,4] = ", L(0,4));

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    ASSERT_EQ(k, e.size());

    //DEBUG(" Static algorithm on L: ");
    Lanczos<CSRMatrix, double> l(L, k, a, skip);
    l.run();
    e = l.getkEigenvalues();
    ASSERT_EQ(k, e.size());
    for (int i = 0; i < e.size(); i++)
        INFO(e[i]);
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    ASSERT_TRUE(l.checkEigenvectors());
    for (int i =0; i < evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        //INFO(" eigenvector : ", i );
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            // INFO(" (comp vs v) : ", evectors[i][j] , ",", v[i][j] );
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-6);
        }
    }    
    
    //DEBUG(" Dynamic algorithm on L: ");
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
    dyn_l.run();
    std::vector<double> dyn_e = dyn_l.getkEigenvalues();
    ASSERT_EQ(k, dyn_e.size());
    
    double err1= 0;
    for(count i=0; i<k; i++) {
        double x = dyn_e[i]-e[i];
        if (x > err1)
            err1 = x;
    }
    if (err1)
        WARN(" Static and dynamic should produce the same result on a static graph.");

    // INFO("dynamic eigen:");
    // for (int i = 0; i < dyn_e.size(); i++)
    //     INFO(dyn_e[i]);
    
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(eigens[i], dyn_e[i], 1e-5);
        EXPECT_NEAR(e[i], dyn_e[i], 1e-5);
    }

    std::vector<Vector> dyn_evectors(k,Vector(n,0.));
    dyn_l.computekEigenvectors();
    dyn_evectors = dyn_l.getkEigenvectors();
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            //INFO(" (comp vs v) : ", dyn_evectors[i][j] , ",", v[i][j] );
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(v[i][j]), 1e-5);
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(evectors[i][j]), 1e-5);
        }
    }    
    // /**** Random edge insertion ****/
    // DEBUG(" Inserting edges ... ");
    // std::vector<GraphEvent> batch;
    // count nInsertions = 3, i = 0;
    // while (i < nInsertions) {
    //     node v1 = GraphTools::randomNode(G);
    //     node v2 = GraphTools::randomNode(G);
    //     if (v1 != v2 && !G.hasEdge(v1, v2)) {
    //         G.addEdge(v1, v2);
    //         L.setValue(v1, v2, -1.0);
    //         L.setValue(v2, v1, -1.0);
    //         L.setValue(v2, v2, L(v2,v2) + 1);
    //         L.setValue(v1, v1, L(v1,v1) + 1);
    //         batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
    //         i++;
    //     }
    // }

    DEBUG(" Inserting edges (1,3) and (3,4) to create updated L'.");
    DEBUG("******************************************** ");
    DEBUG("******************************************** ");
    DEBUG("************  Original Graph  ************** ");
    DEBUG("******************************************** ");
    DEBUG("******************************************** ");
    std::cout << "[ ";
    for (int i = 0; i < eigens.size(); i++) {
        std::cout << eigens[i] << " ";
    }
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< v.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j< n; j++) {
            std::cout << v[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");



    
    std::vector<GraphEvent> batch;
    node v1 = 1;
    node v2 = 3;
    G.addEdge(v1, v2);
    L.setValue(v1, v2, -1.0);
    L.setValue(v2, v1, -1.0);
    L.setValue(v2, v2, L(v2,v2) + 1);
    L.setValue(v1, v1, L(v1,v1) + 1);
    batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));

    v1 = 3;
    v2 = 4;
    G.addEdge(v1, v2);
    L.setValue(v1, v2, -1.0);
    L.setValue(v2, v1, -1.0);
    L.setValue(v2, v2, L(v2,v2) + 1);
    L.setValue(v1, v1, L(v1,v1) + 1);
    batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));


    

    std::vector<double> dyn_eigens = {5.30277563, 4.860805853111, 3.25410168836505, 1.69722436226801, 0.885092458523244, -8.06487158886027e-16 };
    const std::vector<std::vector<double>> dyn_v = {
                                                    { 0.173523521686949, 0.0525386949402678, -0.746632781688064, 0.573109260001114, 0.173523521686949, -0.226062216627217},
                                                    {-0.0934131751839204, -0.360650133507829, 0.360650133507829, 0.671099879356886, -0.527486316514192, -0.0502003876587728},
                                                    {-0.147425200856250, -0.332311394157631, 0.332311394157631, 0.0844408863183978, 0.643166216534716, -0.580181901996864},
                                                    {0.230700933580019, -0.761953423030114,-0.160850311289962,-0.0698506222900568,0.230700933580019, 0.531252489450094},
                                                    {0.848256912058224 ,-0.0974711163052746, 0.0974711163052753, -0.206142398950184,-0.241029881547304,-0.401084631560734},
                                                    {-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863}};

                
        
    DEBUG("******************************************** ");
    DEBUG("******************************************** ");
    DEBUG("*************  Updated Graph  ************** ");
    DEBUG("******************************************** ");
    DEBUG("******************************************** ");
    std::cout << "[ ";
    for (int i = 0; i < dyn_eigens.size(); i++) {
        std::cout << dyn_eigens[i] << " ";
    }
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< dyn_v.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << dyn_v[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");

    
    DEBUG("********** Running static on L' ************");
    l.setup(L);
    l.run();
    e = l.getkEigenvalues();
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    std::cout << "[ ";
    for (int i = 0; i < e.size(); i++)
        std::cout << e[i] << " ";
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< evectors.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << evectors[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");

    
    
    // for (int i = 0; i < e.size(); i++) {
    //     EXPECT_NEAR(dyn_eigens[i], e[i], 1e-5);
    // }
    
    DEBUG("********** Running dynamic on L' *********** ");
    // for dynamic algo the update performed on eigenvalues and eigenvectors

    dyn_l.updateBatch(batch);
    dyn_e = dyn_l.getkEigenvalues();
    dyn_evectors = dyn_l.getkEigenvectors();

    //dyn_e = dyn_l.getkEigenvalues();
    //dyn_evectors = dyn_l.getkEigenvectors();


    err1 = 0;
    for(count i=0; i<k; i++) {
      double x = dyn_e[i]-e[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG(" Static vs Dynamic error = ", err1);
    if (err1)
        WARN(" Static and dynamic should produce the same result on a static graph.");


    std::cout << "[ ";
    for (int i = 0; i < dyn_e.size(); i++)
        std::cout << dyn_e[i] << " ";
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< dyn_evectors.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << dyn_evectors[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");


    
    // DEBUG(" Checking compared to matlab result: ");
    // /* check if result from dymanic is equal to the matlab resulting */
    // // for (int i = 0; i < e.size(); i++) 
    // //     EXPECT_NEAR(e[i], dyn_e[i], 1e-5);

    // assert(e.size() == dyn_e.size());
    // INFO("real \t static \t dynamic ");
    // for (int i = 0; i < e.size(); i++) {
    //     std::cout << "\t" << dyn_eigens[i] << "\t" << e[i] << "\t" <<  dyn_e[i] << std::endl;
    // }

    
    // assert(dyn_evectors.size() == evectors.size());
    // for (int i =0; i< dyn_evectors.size(); i++) {
    //     ASSERT_EQ(n, dyn_evectors[i].getDimension());
    //     ASSERT_EQ(n, evectors[i].getDimension());
    //     INFO("  \t\t Eigenvector ", i);
    //     INFO("real \t static \t dynamic ");
    //     for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
    //         std::cout << "\t" << dyn_v[i][j] << "\t" << evectors[i][j] << "\t" <<  dyn_evectors[i][j] << std::endl;
    //     }
    // }

    // ASSERT_TRUE(l.checkEigenvectors());
    // ASSERT_TRUE(dyn_l.checkEigenvectors());
    // for (int i =0; i< dyn_evectors.size(); i++) {
    //     ASSERT_EQ(n, dyn_evectors[i].getDimension());
    //     for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
    //         //INFO(" (comp vs v) : ", dyn_evectors[i][j] , ",", evectros[i][j] );
    //         EXPECT_NEAR(dyn_evectors[i][j], evectors[i][j], 1e-5);
    //     }
    // }    



    std::vector<GraphEvent> batch2;
    v1 = 0;
    v2 = 1;
    G.addEdge(v1, v2);
    L.setValue(v1, v2, -1.0);
    L.setValue(v2, v1, -1.0);
    L.setValue(v2, v2, L(v2,v2) + 1);
    L.setValue(v1, v1, L(v1,v1) + 1);
    batch2.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));

    DEBUG("************** ONE MORE RUN ****************");
    DEBUG("********** Running static on L' ************");
    
    l.setup(L);
    l.run();
    e = l.getkEigenvalues();
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    std::cout << "[ ";
    for (int i = 0; i < e.size(); i++)
        std::cout << e[i] << " ";
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< evectors.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << evectors[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");

    
    
    // for (int i = 0; i < e.size(); i++) {
    //     EXPECT_NEAR(dyn_eigens[i], e[i], 1e-5);
    // }
    
    DEBUG("********** Running dynamic on L' *********** ");
    // for dynamic algo the update performed on eigenvalues and eigenvectors

    dyn_l.updateBatch(batch2);
    dyn_e = dyn_l.getkEigenvalues();
    dyn_evectors = dyn_l.getkEigenvectors();

    //dyn_e = dyn_l.getkEigenvalues();
    //dyn_evectors = dyn_l.getkEigenvectors();


    err1 = 0;
    for(count i=0; i<k; i++) {
      double x = dyn_e[i]-e[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG(" Static vs Dynamic error = ", err1);
    if (err1)
        WARN(" Static and dynamic should produce the same result on a static graph.");


    std::cout << "[ ";
    for (int i = 0; i < dyn_e.size(); i++)
        std::cout << dyn_e[i] << " ";
    std::cout << "]\n";
    DEBUG("******************************************** ");
    for (int i =0; i< dyn_evectors.size(); i++) {
        std::cout << "************ Eigenvector # " <<  i << "[ ";
        for (int j = 0; j < n; j++) {
            std::cout << dyn_evectors[i][j] << " ";
        }
        std::cout << "]\n";
    }    
    DEBUG("******************************************** ");


    
    
    
    DEBUG("Finished test.");
}



    
  


  
} /* namespace NetworKit */
