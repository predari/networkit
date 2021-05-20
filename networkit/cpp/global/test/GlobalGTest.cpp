/*
 * GlobalGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <networkit/global/ClusteringCoefficient.hpp>
#include <networkit/global/TriangleCounting.hpp>
#include <networkit/global/DynTriangleCounting.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/algebraic/algorithms/AlgebraicTriangleCounting.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>

#include <networkit/generators/ErdosRenyiGenerator.hpp>

namespace NetworKit {

class GlobalGTest: public testing::Test {};

TEST_F(GlobalGTest, testClusteringCoefficient) {

    ErdosRenyiGenerator graphGen(10, 1.0);
    Graph G = graphGen.generate();

    ClusteringCoefficient clusteringCoefficient;
    double cc = clusteringCoefficient.avgLocal(G);

    EXPECT_EQ(1.0, cc);
}



TEST_F(GlobalGTest, testGlobalClusteringCoefficient) {
    Graph G(6);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(2, 5);
    G.addEdge(3, 5);

    double ccg = ClusteringCoefficient::exactGlobal(G);
    EXPECT_NEAR(ccg, 18.0 / 34.0, 1e-9);
}



TEST_F(GlobalGTest, testToyGraphTwo) {
    Graph graph(5);

    graph.addEdge(0,1);
    graph.addEdge(0,2);
    graph.addEdge(1,2);
    graph.addEdge(1,3);
    graph.addEdge(3,4);
    graph.addEdge(2,4);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    EXPECT_TRUE(atc.hasFinished());
    std::vector<count> nodeScores = atc.getScores();

    TriangleCounting tc(graph);
    //tc.run_seq();
    tc.run_par();
    EXPECT_TRUE(tc.hasFinished());
    std::vector<count> nodeTriangles = tc.getTriangleScores();
    INFO("Printing triangle scores:");
    for (count i = 0; i < nodeTriangles.size(); i ++) {
            INFO(" ", nodeTriangles[i]);
    }

    count triangles = tc.getTriangleCount();
    INFO("Printing triangle count = ", triangles);

    EXPECT_EQ(1u, nodeScores[0]) << "wrong algebraic triangle count";
    EXPECT_EQ(1u, nodeScores[1]) << "wrong algebraic triangle count";
    EXPECT_EQ(1u, nodeScores[2]) << "wrong algebraic triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong algebraic triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong algebraic triangle count";

    EXPECT_EQ(1u, nodeTriangles[0]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeTriangles[1]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeTriangles[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeTriangles[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeTriangles[4]) << "wrong triangle count";


    // graph.addEdge(2,3);
    // atc = AlgebraicTriangleCounting<CSRMatrix>(graph);
    // atc.run();

    // EXPECT_EQ(1u, atc.score(0)) << "wrong triangle count";
    // EXPECT_EQ(2u, atc.score(1)) << "wrong triangle count";
    // EXPECT_EQ(3u, atc.score(2)) << "wrong triangle count";
    // EXPECT_EQ(2u, atc.score(3)) << "wrong triangle count";
    // EXPECT_EQ(1u, atc.score(4)) << "wrong triangle count";
}

        

TEST_F(GlobalGTest, testTriangleCounting) {
        // METISGraphReader reader;
        // Graph G = reader.read("../input/PGPgiantcompo.graph");
        SNAPGraphReader reader;
        Graph G = reader.read("../input/wiki-Vote.txt");
        
        // Graph G = reader.read("../input/celegans_metabolic.graph");
        // count n = G.upperNodeIdBound();

        Aux::Timer timer;
        
        TriangleCounting tc(G);
        INFO("Running TriangleCounting. ");
        timer.start();
        tc.run_seq();
        //tc.run_par();
        timer.stop();
        
        count triangles = tc.getTriangleCount()/6;
        printf("\n%s <%d> %f\n", __FUNCTION__, __LINE__, timer.elapsedMicroseconds() / 1e6);
        printf("Triangles in original graph = %ld\n", triangles);
        EXPECT_EQ(608389, triangles); // for wiki-Vote.txt
        
        
}


TEST_F(GlobalGTest, testToyDynTriangleCounting) {

/* Graph:
7 - 2 
 \ / (\) 
   1(-) 3
  /\(\) (/)
6   5 - 4
*/
    int n = 7;
    Graph G(n);

    G.addEdge(0, 1);
    G.addEdge(0, 4);
    G.addEdge(0, 5);
    G.addEdge(0, 6);
    G.addEdge(1, 6);
    G.addEdge(3, 4);
    DEBUG(" ** Original G (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");

    DynTriangleCounting dyntc(G);
    DEBUG(" ** Counting triangles in original graph (static algo).");
    dyntc.run();
    count triangles = dyntc.getTriangleCount();
    INFO("** triangles = ", triangles, ".");
    EXPECT_EQ(1, triangles); // for wiki-Vote.txt
    
    // make one batch with edge in paranthesis (see figure)
    DEBUG(" ** Creating edge addition batch. ");
    std::vector<GraphEvent> addition;
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 2));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 3));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 3, 2));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 2, 1));
    G.addEdge(0,2);
    G.addEdge(0,3);
    G.addEdge(3,2);
    G.addEdge(2,1);
    DEBUG(" ** Updated-Addition G (n, m) = (" ,G.numberOfNodes() , ", " , G.numberOfEdges() , ")");    
    DEBUG(" ** Counting triangles in new graph (static algo). ");
    dyntc.reset(G, true);
    dyntc.run();
    count new_triangles = dyntc.getTriangleCount();
    INFO("** New triangles = ", new_triangles, ".");
    EXPECT_GE(new_triangles, triangles);

    G.removeEdge(0,2);
    G.removeEdge(0,3);
    G.removeEdge(3,2);
    G.removeEdge(2,1);
    DEBUG(" ** Updated-Removal G (n, m) = (" ,G.numberOfNodes() , ", " , G.numberOfEdges() , ")");

    dyntc.reset(G, true);
    assert(dyntc.Insertion());

    dyntc.edgeInsertion(addition);
    assert(dyntc.checkSorted());
    
    dyntc.updateBatch(addition);
    count update_triangles = dyntc.getTriangleCount();
    printf(" ** Triangles created by batch update = %ld ", update_triangles);
    printf(" ... should be equal %ld ==", (new_triangles - triangles)); 
    printf("(new %ld - old %ld)\n", new_triangles, triangles);
    EXPECT_EQ(update_triangles, (new_triangles - triangles)); // for wiki-Vote.txt

        
}        
        



TEST_F(GlobalGTest, testDynTriangleCounting) {

        SNAPGraphReader reader;
        Graph G = reader.read("../input/wiki-Vote.txt");
        Aux::Timer timer;
        DEBUG(" ** Original graph:  G (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");

        DynTriangleCounting dyntc(G);
        DEBUG(" ** Counting triangles in original graph (static algo).");
        timer.start();
        dyntc.run();
        timer.stop();
        
        count triangles = dyntc.getTriangleCount();
        INFO("** triangles = ", triangles, " in ", timer.elapsedMicroseconds() / 1e6, " secs.");
        EXPECT_EQ(608389, triangles); // for wiki-Vote.txt
        
	// Make multiple batches of removed edges
        DEBUG(" ** Creating edge removal batches. ");
        count numBatches = 5;
        count numEdges = 100;
        //count dupEdges = 10;
        //count totalEdges = numEdges + dupEdges;
        std::vector<std::vector<GraphEvent>> deletions(numBatches);
        
        for (count i = 0; i < numBatches; i++) {
                for (count j = 0; j < numEdges; j++) {
                        
                        node u = G.upperNodeIdBound();
                        node v = G.upperNodeIdBound();
                        do {
                                u = GraphTools::randomNode(G);
                                v = GraphTools::randomNode(G);
                        } while (!G.hasEdge(u, v) || u == v);
                        //DEBUG(" ** Edge (", u, ", ", v, ") # ", j , " to be removed. ");       
                        GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
                        deletions[i].push_back(edgeDeletion);
                        G.removeEdge(u, v);
                }
                // // Add duplicates from beginning of batch
                // for(unsigned j = 0; j < dupEdges; ++j) {
                //         assert(j < deletions.size());
                //         auto dupEvent = deletions[i][j];
                //         deletions[i].push_back(GraphEvent::EDGE_REMOVAL, dupEvent.u, dupEvent.v);
                // }            
        }
        // updated graph has been created.
        DEBUG(" ** Graph after removing edges:  G new (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
        
        DEBUG(" ** Counting triangles in new graph (static algo). ");
        timer.start();
        dyntc.reset(G, false);
        dyntc.run();
        timer.stop();       
        count new_triangles = dyntc.getTriangleCount();
        INFO("** New triangles = ", new_triangles, " in ", timer.elapsedMicroseconds() / 1e6, " secs.");
        EXPECT_LE(new_triangles, triangles);
        // TODO: make sure sorting is done only once. Check member functions run() and updateBatch() 
        G.sortEdges();

        dyntc.reset(G, true);
        if (!dyntc.Insertion())
                WARN("** Update test is made on insertion mode!!");
        
        count total_update_triangles = 0;
        count total_update_triangles_ = 0;
        for(unsigned i = 0; i < numBatches; ++i) {
                assert(i < deletions.size());
                timer.start();
                // check insertion or deletion in general
                // dyntc.edgeInsertionSorted(deletions[i]);
                dyntc.edgeInsertion(deletions[i]);
                if(!dyntc.checkSorted())
                        WARN("** Graph after insertion is not sorted!!");
                timer.stop();
                INFO("** Time for inserting (sorted) batch into (sorted) graph in ", timer.elapsedMicroseconds() / 1e6, " secs.");
                timer.start();
                dyntc.updateBatch(deletions[i]);
                timer.stop();       
                INFO("** Time for updateBatch # ", i, "in ", timer.elapsedMicroseconds() / 1e6, " secs.");
                total_update_triangles_ =+ dyntc.getTriangleCount();
        }

        total_update_triangles = dyntc.getTriangleCount();
        printf(" ** Triangles created by batch update = %ld ", total_update_triangles);
        
	// Let's compare
	printf(" == %ld  == %ld ", total_update_triangles_, triangles);
	printf("old %ld, new %ld\n", triangles, new_triangles);
        EXPECT_EQ(total_update_triangles, total_update_triangles_); // for wiki-Vote.txt
        EXPECT_EQ(total_update_triangles, triangles - new_triangles); // for wiki-Vote.txt

        
}        


        
} /* namespace NetworKit */
