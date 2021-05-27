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
        SNAPGraphReader reader;
        Graph G = reader.read("../input/wiki-Vote.txt");
        
        Aux::Timer timer;
        
        TriangleCounting tc(G);
        INFO("Running TriangleCounting. ");
        timer.start();
        tc.run_seq();
        timer.stop();
        
        count triangles = tc.getTriangleCount()/6;
        printf("\n%s <%d> %f\n", __FUNCTION__, __LINE__, timer.elapsedMicroseconds() / 1e6);
        printf("Triangles in original graph = %ld\n", triangles);
        EXPECT_EQ(608389, triangles); // for wiki-Vote.txt
        
        
}


TEST_F(GlobalGTest, testToyDynTriangleCountingI) {

/* Graph: 
7 - 2 
 \ / (\) 
   1(-) 3
  /\(\) (/)
6   5 - 4
* edges in () correspond to updates.
*/
    int n = 7;
    Graph G(n);
    G.addEdge(0, 1);
    G.addEdge(0, 4);
    G.addEdge(0, 5);
    G.addEdge(0, 6);
    G.addEdge(1, 6);
    G.addEdge(3, 4);
    INFO(" ** STATIC ALGO on G (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
    DynTriangleCounting dyntc(G);
    dyntc.run();
    double triangles = dyntc.getTriangleCount();
    EXPECT_EQ(1, triangles);
    INFO(" ** ***** Triangles = ", triangles, ".");

    

    std::vector<GraphEvent> addition;
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 2));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 3));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 3, 2));
    addition.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 2, 1));
    
    INFO(" ** DYNAMIC ALGO on updated G_I = {G}U{G'}.");
    //dyntc.edgeInsertion(addition);
    assert(dyntc.checkSorted());
    dyntc.updateBatch(addition);
    double final_triangles = dyntc.getTriangleCount();
    EXPECT_EQ(dyntc.getNewTriangles(), (final_triangles - triangles));
    EXPECT_EQ(4, final_triangles);
    INFO(" ** ***** Triangles  = ", final_triangles, " = (", triangles, " + ",
         dyntc.getNewTriangles(), ")." );
    
    /* checking triangle counting score calculation*/
    std::vector<double> t_score = dyntc.getTriangleScores();
    assert(t_score.size() == G.numberOfNodes());
    double t_t = 0;
    for (node u = 0; u < G.upperNodeIdBound(); u++)
            t_t += t_score[u];
    EXPECT_EQ(t_t/3, final_triangles);

    INFO(" ** ***** TEST ");
    INFO(" ** ***** G = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
    INFO(" ** STATIC ALGO on G_I. ");
    dyntc.run();
    double static_triangles = dyntc.getTriangleCount();
    EXPECT_EQ(static_triangles, final_triangles);
    EXPECT_EQ(0, dyntc.getNewTriangles()); 
    INFO("** ***** Triangles = ", static_triangles, ".");
    
}        
        



TEST_F(GlobalGTest, testToyDynTriangleCountingD) {

/* Graph: 
6 - 1 
 \ / (\) 
   0(-) 2
  /\(\) (/)
5   4 - 3
* edges in () correspond to updates (deletion).
*/
    int n = 7;
    Graph G(n);

    G.addEdge(0, 1);
    G.addEdge(0, 4);
    G.addEdge(0, 5);
    G.addEdge(0, 6);
    G.addEdge(1, 6);
    
    G.addEdge(3, 4); 
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(3, 2);
    G.addEdge(2, 1);

    

    DynTriangleCounting dyntc(G);
    INFO(" ** STATIC ALGO on G (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
    dyntc.run();
    double triangles = dyntc.getTriangleCount();
    EXPECT_EQ(4, triangles);
    INFO(" ** ***** Triangles = ", triangles, ".");
    
    std::vector<double> t_score = dyntc.getTriangleScores();
    assert(t_score.size() == G.numberOfNodes());
    double t_t = 0;
    for (node u = 0; u < G.upperNodeIdBound(); u++)
            t_t += t_score[u];
    EXPECT_EQ(t_t/3, triangles);
   
    
    std::vector<GraphEvent> deletion;
    deletion.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, 0, 2));
    deletion.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, 0, 3));
    deletion.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, 3, 2));
    deletion.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, 2, 1));
    

    //dyntc.edgeDeletion(deletion);
    assert(dyntc.checkSorted());
    INFO(" ** DYNAMIC ALGO on updated G_D = {G}\\cap{G'}.");
    dyntc.updateBatch(deletion);
    double final_triangles = dyntc.getTriangleCount();
    INFO(" ** ***** Triangles  = ", final_triangles, " = (", triangles, " - ",
         dyntc.getNewTriangles(), ")." );
    EXPECT_EQ(1, final_triangles);
    EXPECT_EQ(3, dyntc.getNewTriangles()); 

    t_score = dyntc.getTriangleScores();
    assert(t_score.size() == G.numberOfNodes());
    t_t = 0;
    for (node u = 0; u < G.upperNodeIdBound(); u++) 
            t_t += t_score[u];
    EXPECT_EQ(t_t/3, final_triangles);
    
    /* run static algorithm on the new graph */
    INFO(" ** STATIC ALGO on G_D.");
    dyntc.run();
    double static_triangles = dyntc.getTriangleCount();
    EXPECT_EQ(static_triangles, final_triangles);
    EXPECT_EQ(0, dyntc.getNewTriangles());
    INFO("** ***** Triangles = ", static_triangles, ".");
    
}        




        

TEST_F(GlobalGTest, testDynTriangleCounting) {

        SNAPGraphReader reader;
        Graph G = reader.read("../input/wiki-Vote.txt");
        Aux::Timer timer;

        INFO(" ** STATIC ALGO on final G (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
        DynTriangleCounting dyntc(G);
        timer.start();
        dyntc.run();
        timer.stop();
        
        double triangles = dyntc.getTriangleCount();
        EXPECT_EQ(608389, triangles); // for wiki-Vote.txt
        INFO(" ** ***** Triangles  = ", triangles, " in ", timer.elapsedMicroseconds() / 1e6, " secs.");
        double static_time = timer.elapsedMicroseconds() / 1e6;

        
	// Make multiple batches of removed edges
        count numBatches = 5;
        count numEdges = 100;
        count dupEdges = 10;
        //count totalEdges = numEdges + dupEdges;
        std::vector<std::vector<GraphEvent>> addition(numBatches);
        
        for (count i = 0; i < numBatches; i++) {
                
                for (count j = 0; j < numEdges; j++) {                        
                        node u = G.upperNodeIdBound();
                        node v = G.upperNodeIdBound();
                        do {
                                u = GraphTools::randomNode(G);
                                v = GraphTools::randomNode(G);
                        } while (!G.hasEdge(u, v) || u == v);
                        GraphEvent edgeI(GraphEvent::EDGE_ADDITION, u, v);
                        addition[i].push_back(edgeI);
                        G.removeEdge(u, v);
                }
                // // Add duplicates from beginning of batch
                // for(unsigned j = 0; j < dupEdges; ++j) {
                //         assert(j < numEdges);
                //         auto dupEvent = addition[i][j];
                //         addition[i].push_back(dupEvent);
                // }            
        }
        
        INFO(" ** STATIC ALGO on start G (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
        timer.start();
        dyntc.run();
        timer.stop();       
        double static_triangles = dyntc.getTriangleCount();
        INFO(" ** ***** Triangles = ", static_triangles, " in ", timer.elapsedMicroseconds() / 1e6, " secs.");
        double starting_time = timer.elapsedMicroseconds() / 1e6;

        EXPECT_LE(static_triangles, triangles);
        EXPECT_EQ(0, dyntc.getNewTriangles()); 
        // TODO: make sure sorting is done only once.
        // TODO: input should be sorted for edgeInsertionSorted.
        //G.sortEdges();
        INFO(" ** DYNAMIC ALGO on updated G_I = {G}U{G'}.");
        double total_update_triangles = 0;
        double total_update_time = 0.0;
        for(unsigned i = 0; i < numBatches; ++i) {
                assert(i < addition.size());
                timer.start();
                //dyntc.edgeInsertion(addition[i]);
                timer.stop();
                if(!dyntc.checkSorted())
                        WARN(" ** ***** GRAPH IS NOT SORTED!");
                INFO(" ** ***** Time to insert batch[",i,"] =  ", timer.elapsedMicroseconds() / 1e6, " secs.");
                total_update_time += timer.elapsedMicroseconds() / 1e6;
                timer.start();
                dyntc.updateBatch(addition[i]);
                timer.stop();       
                INFO(" ** ***** Time to update batch[",i,"] = ", timer.elapsedMicroseconds() / 1e6, " secs.");
                total_update_triangles += dyntc.getNewTriangles();
                total_update_time += timer.elapsedMicroseconds() / 1e6;
                
        }


        
        double final_triangles = dyntc.getTriangleCount();
        EXPECT_EQ(final_triangles, triangles);        
        EXPECT_EQ(total_update_triangles, triangles - static_triangles);
        INFO(" ** ***** Triangles  = ", final_triangles, " = (", static_triangles, " + ",
             total_update_triangles, ")." );
        
        INFO(" ** STATIC_TIME   =  ", static_time );
        INFO(" ** DYNAMIC_TIME  =  ", total_update_time + starting_time );
        

        std::vector<double> t_score = dyntc.getTriangleScores();
        assert(t_score.size() == G.numberOfNodes());
        double t_t = 0;

        for (node u = 0; u < G.upperNodeIdBound(); u++) {
                t_t += t_score[u];
        }
        EXPECT_EQ(t_t/3, final_triangles);
        
        
}        


        
} /* namespace NetworKit */
