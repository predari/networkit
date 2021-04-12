#include <gtest/gtest.h>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/DynLanczos.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class DynLanczosGTest : public testing::Test {};


TEST(DynLanczosGTest, runDynVsStOnStSmallGraph) {

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

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);
    ASSERT_EQ(k, e.size());

    // runing static algorithm
    Lanczos<CSRMatrix, double> l(L, k, a, skip);
    l.run();
    e = l.getkEigenvalues();
    ASSERT_LE(e.size(), k);
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(eigens[i], e[i], 1e-5);
    }

    std::vector<Vector> evectors(k,Vector(n,0.));
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    ASSERT_TRUE(l.checkEigenvectors());
    for (int i =0; i < evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(v[i][j]), 1e-6);
        }
    }    
    
    // runing dynamic algorithm
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
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
        }
    }
}



TEST(DynLanczosGTest, runDynVsStOnSt) {
    METISGraphReader reader;
    Graph G = reader.read("../input/PGPgiantcompo.graph");
    //Graph G = reader.read("input/celegans_metabolic.graph");
    count n = G.upperNodeIdBound();

    auto L = CSRMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;
    std::vector<double> e(k,0.);

    // runing static algorithm
    Lanczos<CSRMatrix, double> l(L, k, a, skip);
    l.run();
    e = l.getkEigenvalues();
    ASSERT_LE(e.size(), k);

    std::vector<Vector> evectors(k,Vector(n,0.));
    l.computekEigenvectors();
    evectors = l.getkEigenvectors();
    ASSERT_TRUE(l.checkEigenvectors());

    
    // runing dynamic algorithm
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
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
        EXPECT_NEAR(e[i], dyn_e[i], 1e-5);
    }

    std::vector<Vector> dyn_evectors(k,Vector(n,0.));
    dyn_l.computekEigenvectors();
    dyn_evectors = dyn_l.getkEigenvectors();
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(evectors[i][j]), 1e-5);
        }
    }

    

}

    

TEST(DynLanczosGTest, DynVsStOnDynSmallGraph) {

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

    // original graph
    auto L = CSRMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;

    Lanczos<CSRMatrix, double> l(L, k, a, skip);
    l.run();
    l.computekEigenvectors();
    
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
    dyn_l.run();
    dyn_l.computekEigenvectors();


    DEBUG(" Inserting edges (1,3) and (3,4) to create updated L'.");
    
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

    // computed from matlab
    std::vector<double> dyn_eigens = {5.30277563, 4.860805853111, 3.25410168836505, 1.69722436226801, 0.885092458523244, -8.06487158886027e-16 };
    const std::vector<std::vector<double>> dyn_v = {
                                                    { 0.173523521686949, 0.0525386949402678, -0.746632781688064, 0.573109260001114, 0.173523521686949, -0.226062216627217},
                                                    {-0.0934131751839204, -0.360650133507829, 0.360650133507829, 0.671099879356886, -0.527486316514192, -0.0502003876587728},
                                                    {-0.147425200856250, -0.332311394157631, 0.332311394157631, 0.0844408863183978, 0.643166216534716, -0.580181901996864},
                                                    {0.230700933580019, -0.761953423030114,-0.160850311289962,-0.0698506222900568,0.230700933580019, 0.531252489450094},
                                                    {0.848256912058224 ,-0.0974711163052746, 0.0974711163052753, -0.206142398950184,-0.241029881547304,-0.401084631560734},
                                                    {-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863,-0.408248290463863}};

                
        
    
    DEBUG("******************************************** ");
    DEBUG("*********** Updated Graph Eigen ************ ");
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

    
    // setup again, with new input graph 
    l.setup(L);
    l.run();
    std::vector<double> e = l.getkEigenvalues();
    l.computekEigenvectors();
    std::vector<Vector> evectors = l.getkEigenvectors();
    ASSERT_LE(evectors.size(), k);
    // std::cout << "[ ";
    // for (int i = 0; i < e.size(); i++)
    //     std::cout << e[i] << " ";
    // std::cout << "]\n";
    // DEBUG("******************************************** ");
    // for (int i =0; i< evectors.size(); i++) {
    //     std::cout << "************ Eigenvector # " <<  i << "[ ";
    //     for (int j = 0; j < n; j++) {
    //         std::cout << evectors[i][j] << " ";
    //     }
    //     std::cout << "]\n";
    // }    
    // DEBUG("******************************************** ");
    
    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(dyn_eigens[i], e[i], 1e-5);
    }
    ASSERT_TRUE(l.checkEigenvectors());
    for (int i =0; i< evectors.size(); i++) {
        ASSERT_EQ(n, evectors[i].getDimension());
        for (int j = 0; j< evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(evectors[i][j]), fabs(dyn_v[i][j]), 1e-5);
        }
    }

    // running dynamic update on eigenvalues and eigenvectors simmultaneously
    dyn_l.updateBatch(batch);
    std::vector<double> dyn_e = dyn_l.getAllEigenvalues();
    std::vector<Vector> dyn_evectors = dyn_l.getBasis();

    double err1= 0;
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

    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(dyn_eigens[i], dyn_e[i], 1e-5);
    }
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(dyn_v[i][j]), 1e-5);

        }
    }    

}

    

TEST(DynLanczosGTest, DynVsStOnDyn) {

    METISGraphReader reader;
    Graph G = reader.read("../input/PGPgiantcompo.graph");
    //Graph G = reader.read("input/celegans_metabolic.graph");
    count n = G.upperNodeIdBound();

    auto L = CSRMatrix::laplacianMatrix(G);

    int k = 2;
    int a = 2;
    int skip = 16;

    Lanczos<CSRMatrix, double> l(L, k, a, skip);
    l.run();
    l.computekEigenvectors();
    
    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
    dyn_l.run();
    dyn_l.computekEigenvectors();


    DEBUG(" Inserting random edges to create updated L'.");


    /**** Random edge insertion ****/
    std::vector<GraphEvent> batch;
    count nInsertions = 6, i = 0;
    while (i < nInsertions) {
        node v1 = GraphTools::randomNode(G);
        node v2 = GraphTools::randomNode(G);
        if (v1 != v2 && !G.hasEdge(v1, v2)) {
            G.addEdge(v1, v2);
            L.setValue(v1, v2, -1.0);
            L.setValue(v2, v1, -1.0);
            L.setValue(v2, v2, L(v2,v2) + 1);
            L.setValue(v1, v1, L(v1,v1) + 1);
            batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
            i++;
        }
    }
    
    // setup again, with new input graph 
    l.setup(L);
    l.run();
    std::vector<double> e = l.getkEigenvalues();
    l.computekEigenvectors();
    std::vector<Vector> evectors = l.getkEigenvectors();
    ASSERT_LE(evectors.size(), k); 
    //ASSERT_TRUE(l.checkEigenvectors());
  

    // running dynamic update on eigenvalues and eigenvectors simmultaneously
    dyn_l.updateBatch(batch);
    std::vector<double> dyn_e = dyn_l.getAllEigenvalues();
    std::vector<Vector> dyn_evectors = dyn_l.getBasis();
    
    double err1= 0;
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

    for (int i = 0; i < e.size(); i++) {
        EXPECT_NEAR(e[i], dyn_e[i], 1e-5);
    }
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(evectors[i][j]), 1e-5);

        }
    }    

}


TEST(DynLanczosGTest, DynBatchSeries) {

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

    int k = 2;
    int a = 2;
    int skip = 16;


    DynLanczos<CSRMatrix,double> dyn_l(L, k, a, skip);
    dyn_l.run();
    std::vector<double> dyn_e = dyn_l.getkEigenvalues();
    ASSERT_LE(dyn_e.size(), k);
        
    for (int i = 0; i < dyn_e.size(); i++) {
        EXPECT_NEAR(eigens[i], dyn_e[i], 1e-5);
    }

    std::vector<Vector> dyn_evectors(k,Vector(n,0.));
    dyn_l.computekEigenvectors();
    dyn_evectors = dyn_l.getkEigenvectors();
    ASSERT_TRUE(dyn_l.checkEigenvectors());
    
    for (int i =0; i< dyn_evectors.size(); i++) {
        ASSERT_EQ(n, dyn_evectors[i].getDimension());
        for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
            EXPECT_NEAR(fabs(dyn_evectors[i][j]), fabs(v[i][j]), 1e-5);
        }
    }    

    DEBUG(" Inserting edges (1,3) and (3,4) to create updated L'.");
    
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

    dyn_l.updateBatch(batch);
    dyn_e = dyn_l.getAllEigenvalues();
    dyn_evectors = dyn_l.getBasis();
    

    // for (int i = 0; i < dyn_e.size(); i++) 
    //     EXPECT_NEAR(dyn_eigens[i], dyn_e[i], 1e-5);

    ASSERT_TRUE(dyn_l.checkEigenvectors());
    // for (int i =0; i< dyn_evectors.size(); i++) {
    //     ASSERT_EQ(n, dyn_evectors[i].getDimension());
    //     for (int j = 0; j< dyn_evectors[i].getDimension(); j++) {
    //         EXPECT_NEAR(dyn_v[i][j], dyn_evectors[i][j], 1e-5);
    //     }
    // }    


    std::vector<GraphEvent> batch2;
    v1 = 0;
    v2 = 1;
    G.addEdge(v1, v2);
    batch2.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
    
    dyn_l.updateBatch(batch2);
    dyn_e = dyn_l.getAllEigenvalues();
    dyn_evectors = dyn_l.getBasis();

    ASSERT_TRUE(dyn_l.checkEigenvectors());

}    


    
  


  
} /* namespace NetworKit */
