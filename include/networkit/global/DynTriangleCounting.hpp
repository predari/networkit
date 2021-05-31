/*
 * DynTriangleCounting.hpp
 *
 *  Created on: May 2021
 *      Author: Predari Maria
 */

#ifndef NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_
#define NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_

#include <omp.h>

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

/**
 * Calculates the number of triangles in a dynamic graph.
 *
 * @ingroup global
 */
class DynTriangleCounting : public Algorithm, public DynAlgorithm {

public:
  /**
   * Constructs the DynTriangleCounting class. This class implements 
   * algorithms for computing the total number of triangles in the graph 
   * (via triangle counts per vertex ). 
   * The implementation is based on the static algorithm in stinger and the 
   * dynamic implementation is taken from Oded Green et al. 
   *
   * @param G The graph.
   */
        
        DynTriangleCounting(Graph &G, bool graphAlreadyUpdated = false)
                : G(&G), TrianglesPerNode(G.upperNodeIdBound(), 0.0),
                  t_t(0.0), o_t(0.0), graphAlreadyUpdated(graphAlreadyUpdated) {
        }
        
        /** 
         * Runs the static algorithm for computing triangle in a graph 
         */
        void run() override;

        
        /**
         * Inserts edges following a dynamic batch update.
         * The method includes sorting of the graph edges after insertion.
         *
         * @param batch A vector of edge modification events.
         */
        void edgeInsertion(const std::vector<GraphEvent> &batch);

        /**
         * Deletes edges following a dynamic batch update.
         * The method includes sorting of the graph edges after deletion.
         *
         * @param batch A vector of edge modification events.
         */
        void edgeDeletion(const std::vector<GraphEvent> &batch);


        /**
         * Updates the total number of triangles in G after a batch of updates.
         *
         * @param batch A vector of edge modification events.
         */
        void updateBatch(const std::vector<GraphEvent> &batch) override;


        /**
         * Updates the total number of triangles in G.
         *
         * @param event The edge modification event.
         */
        
        void update(GraphEvent singleEvent) override {
                std::vector<GraphEvent> events{singleEvent};
                updateBatch(events);
        }

        /**
         * Get total number of triangles
         */
        count getTriangleCount() {
                assureFinished();
                return t_t;
        }
        /**
         * Get new triangles 
         */
        count getNewTriangles() {
                assureFinished();
                return abs(t_t-o_t);
        }
        /**
         * Compute total number of triangles
         */
        count computeTriangleCount() {
                t_t = G->parallelSumForNodes([&](node u) {
                                                     return TrianglesPerNode[u];
                                             });
                
                t_t = t_t/3;
        }
        
        /**
         * Get the triangle count score per node
         */
        
        std::vector<double> getTriangleScores() {
                assureFinished();
                return TrianglesPerNode;
        }

        /**
         * Check if graph is sorted
         */
        
        bool checkSorted(const Graph * G1);

        
private:

        count bsearch_intersection( node u, node v, const std::vector<node> & u_adj, count u_deg);
        
        count sorted_intersection( node u, const std::vector<node> & u_adj, count u_deg,
                                   node v,const std::vector<node> & v_adj, count v_deg, double m);
        // consider call by ref for adj_inG and adj_inB
        count do_merge_insert(const std::vector<node> adj_inG, const count deg_inG,
                              const std::vector<node> adj_inB, const count deg_inB,
                              std::vector<node> & adj_inGI);

        count do_merge_delete(std::vector<node> & adj_inG, const count deg_inG,
                              const std::vector<node> & adj_inB, const count deg_inB);

        void countUpdateTriangles(const Graph &ugraph, bool insertion);

        Graph edgeInsertionSorted(const Graph & ugraph);
        
        Graph edgeDeletionSorted(const Graph & ugraph);

        // pointer to new graph
        Graph *G;
        Graph U;
        // current triangle count per node
        std::vector<double> TrianglesPerNode;
        // current total number of triangles
        double t_t;
        // old total number of triangles
        double o_t;
        bool graphAlreadyUpdated;
        

};

        // from stringer 
        count DynTriangleCounting::bsearch_intersection(node u, node v, const std::vector<node> & u_adj, count u_deg) {
        count t = 0;
        // parallel?
        G->forEdgesOf(v, [&](node, node w) {
                                 if (w != u) {
                                         int64_t first = 0;
                                         int64_t last = u_deg-1;
                                         int64_t middle = (first + last)/2;
                                         
                                         while (first <= last) {
                                                 if (u_adj[middle] < w) {
                                                         first = middle + 1;
                                                 } else if (u_adj[middle] == w) {
                                                         t++;
                                                         break;
                                                 } else {
                                                         last = middle - 1;
                                                 }                 
                                                 middle = (first + last)/2;
                                         }
                                 }     
                         });
        return t;
}
        

        void DynTriangleCounting::run() {
                
                std::fill(TrianglesPerNode.begin(),  TrianglesPerNode.end(), 0.0);
                TrianglesPerNode.resize(G->upperNodeIdBound(), 0.0);
                t_t = 0;
                
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                /* sort edge lists */
                // Graph needs to be sorted. Make the assumption upon calling method.
                if (!checkSorted(NULL)) {
                        G->sortEdges();
                        INFO("** ***** SORTING GRAPH *****");
                }
                else {
                        INFO("** ***** GRAPH ALREADY SORTED UPON RUN *****");
                }
                /* copy edge lists */ 
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });
                
                G->balancedParallelForNodes([&](node s) {                                            
                                                    count s_deg = G->degree(s); 
                                                    count c = 0;
                                                    if (edges[s].empty()) {
                                                            assert(G->isIsolated(s));
                                                            assert(!s_deg);
                                                    }
                                                    G->forEdgesOf(s, [&](node, node t) {
                                                                             if ( t != s ) {
                                                                                     c += bsearch_intersection(s, t, edges[s], s_deg);
                                                                             }
                                                                     });
                                                    TrianglesPerNode[s] = c/2;
                                            });
                hasRun = true;
                computeTriangleCount();
                o_t = t_t;
        }


        void DynTriangleCounting::updateBatch(const std::vector<GraphEvent>& batch) {
                
                o_t = t_t;
                bool insertion = true;
                if (batch[0].type == GraphEvent::EDGE_REMOVAL)
                        insertion = false;
                else if (batch[0].type != GraphEvent::EDGE_ADDITION)
                        throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
                
                Graph ugraph = Graph(G->upperNodeIdBound());
                for(auto e : batch){
                        // check that all events of the batch are of the same type.
                        if (e.type == GraphEvent::EDGE_ADDITION && !insertion) {
                                throw std::runtime_error("All Events of an insert batch should be of type EDGE_ADDITION.");
                        } else if (e.type == GraphEvent::EDGE_REMOVAL && insertion) {
                                throw std::runtime_error("All Events of a remove batch should be of type EDGE_REMOVAL.");
                        }
                        if (e.type != GraphEvent::EDGE_ADDITION && e.type != GraphEvent::EDGE_REMOVAL) {
                                throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
                        }
                        ugraph.addEdge(e.u, e.v);
                }

                ugraph.sortEdges();
                INFO("** ***** SORTING BATCH *****");

                // create local complete copies for G
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });
                                
                // create local complete and reduced copies for update graph
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                std::vector<std::vector<node> > red_edges2(ugraph.upperNodeIdBound());
                std::vector<node> reduced_degree2(ugraph.upperNodeIdBound(), 0);

                ugraph.parallelForNodes([&](node u) {
                                                edges2[u].reserve(ugraph.degree(u));
                                                red_edges2[u].reserve(ugraph.degree(u));
                                                ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                             if (v < u) {
                                                                                     red_edges2[u].emplace_back(v);
                                                                                     reduced_degree2[u]++;
                                                                             }
                                                                             edges2[u].emplace_back(v);
                                                                     });
                                        });

                

                Aux::Timer timer;
                if (graphAlreadyUpdated) INFO("** ***** GRAPH ALREADY UPDATED *****");
                else {
                        assert(checkSorted(NULL));
                        assert(checkSorted(&ugraph));
                        std::vector<edgeweight> empty;
                        
                        timer.start();
                        
                        GraphBuilder result(G->upperNodeIdBound(), false, false);
                        G->balancedParallelForNodes([&](node v) {
                                                            count current_degree = G->degree(v);
                                                            count update_degree = ugraph.degree(v);
                                                            count max_new_degree = current_degree + update_degree;
                                                            // isolated nodes
                                                            if (!current_degree && !update_degree)
                                                                    return;                                                     
                                                            if (insertion) {
                                                                    if (update_degree && current_degree) {
                                                                            // INFO("** ***** MERGING ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < edges[v].size(); i++)
                                                                            //         std::cout << edges[v][i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < edges2[v].size(); i++)
                                                                            //         std::cout << edges2[v][i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            
                                                                            std::vector<node> new_edges(max_new_degree); // reserve ?
                                                                            count new_degree = do_merge_insert(edges[v], current_degree,
                                                                                                               edges2[v], update_degree,
                                                                                                               new_edges);
                                                                            assert(new_edges.size() == new_degree);
                                                                            // INFO("** ***** new_degree = ", new_degree );
                                                                            // INFO("** ***** new_edges.size() = ", new_edges.size() );
                                                                            // INFO("** ***** AFTER MERGING OF ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < new_edges.size(); i++)
                                                                            //         std::cout << new_edges[i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            edges[v] = new_edges;
                                                                            // INFO("** ***** AFTER MOVING OF ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < new_edges.size(); i++)
                                                                            //         std::cout << new_edges[i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                          
                                                                            // TODO: currently only for unweighted graphs
                                                                            result.assignNeighborhood(v, new_edges, empty, false);
                                                                            // store the edgelist of v after merging  
                                                                            
                                                                            // INFO("** ***** AFTER SWAPING OF ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < new_edges.size(); i++)
                                                                            //         std::cout << new_edges[i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            
                                                                            
                                                                    }
                                                                    else if (!update_degree) {
                                                                            // INFO("** ***** DIRECT MERGING G ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < edges[v].size(); i++)
                                                                            //         std::cout << edges[v][i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            result.assignNeighborhood(v, edges[v], empty, false);
                                                                    }
                                                                    else {

                                                                            // INFO("** ***** DIRECT MERGING B ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < edges2[v].size(); i++)
                                                                            //         std::cout << edges2[v][i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            result.assignNeighborhood(v, edges2[v], empty, false);
                                                                            // not move, simply copy
                                                                            edges[v] = edges2[v];
                                                                            
                                                                    }
                                                                    
                                                            }
                                                            else {
                                                                    if (update_degree) {
                                                                            count new_degree = do_merge_delete(edges[v], current_degree,
                                                                                                               edges2[v], update_degree);
                                                                            // INFO("** ***** AFTER MERGING OF ", v);
                                                                            // std::cout << "[";
                                                                            // for (int i = 0; i < edges.size(); i++)
                                                                            //         std::cout << edges[i] << " ";
                                                                            // std::cout << "]" << std::endl;
                                                                            assert(new_degree == edges[v].size());
                                                                    }
                                                                    // just copy previous adj
                                                                    // better move instead of swap
                                                                    result.assignNeighborhood(v, edges[v], empty, false);
                                                                   
                                                            }
                                                    });
                                                            
                        Graph GU = result.toGraph(false,true);
                        U = std::move(GU);
                        // if (insertion)
                        //         INFO(" ** ***** INSERTION RUN SUCCESSFULLY ***** ");
                        // else
                        //         INFO(" ** ***** DELETION RUN SUCCESSFULLY ***** ");

                        G = &U;
                        timer.stop();
                        INFO(" ** ***** Time to insert batch =  ", timer.elapsedMicroseconds() / 1e6, " secs.");  
                }
                // INFO(" ** ***** AFTER GRAPH ASSIGNMENT ***** ");
                // INFO(" ** ***** G = (" , G->numberOfNodes() , ", " , G->numberOfEdges() , ")");
                // G->forNodes([&](node u) {
                //                     G->forEdgesOf(u, [&](node, node v, edgeid) {
                //                                             std::cout << " " << v;
                //                                     });
                //                     std::cout << std::endl;
                //             });
                
                double mtype1, mtype2, mtype3;
                if (insertion) {
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        if (!t_t) return;
                }


                // create reduced copies for updated G
                std::vector<std::vector<node> > red_edges(G->upperNodeIdBound());
                std::vector<node> reduced_degree(G->upperNodeIdBound(), 0);
                G->parallelForNodes([&](node u) {
                                            red_edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     // do not loop over after first v < u
                                                                     if (v < u) {
                                                                             red_edges[u].emplace_back(v);
                                                                             reduced_degree[u]++;        
                                                                     } else return;
                                                             });
                                    });

                // INFO(" ** ***** PRINTIG ***** ");
                // INFO(" ** ***** EDGES ***** ");
                // G->forNodes([&](node u) {
                //                     std::cout << u << " [";
                //                     for (int i = 0; i < edges[u].size(); i++)
                //                             std::cout << edges[u][i] << " ";
                //                     std::cout << "]" << std::endl;
                //             });
                // INFO(" ** ***** RED_EDGES ***** ");
                // G->forNodes([&](node u) {
                //                     std::cout << u << " [";
                //                     for (int i = 0; i < red_edges[u].size(); i++)
                //                             std::cout << red_edges[u][i] << " ";
                //                     std::cout << "]" << std::endl;
                //             });
                



                // INFO(" ** ***** UGRAPH EDGES ***** ");
                // ugraph.forNodes([&](node u) {
                //                         std::cout << u << " [";
                //                         for (int i = 0; i < edges2[u].size(); i++)
                //                                 std::cout << edges2[u][i] << " ";
                //                         std::cout << "]" << std::endl;
                //                 });
                // INFO(" ** ***** UGRAPH RED_EDGES ***** ");
                // ugraph.forNodes([&](node u) {
                //                         std::cout << u << " [";
                //                         for (int i = 0; i < red_edges2[u].size(); i++)
                //                                 std::cout << red_edges2[u][i] << " ";
                //                         std::cout << "]" << std::endl;
                //                 });

                
                ugraph.parallelForEdges([&](node u, node v) {
                                                
                                                // work on updated graph only
                                                sorted_intersection(u,edges[u], G->degree(u), v, edges[v], G->degree(v), mtype1); 
                                                // work on updated and update graph
                                                sorted_intersection(u, red_edges[u], reduced_degree[u],
                                                                    v, edges2[v], ugraph.degree(v), mtype2);
                                                sorted_intersection(v, red_edges[v], reduced_degree[v],
                                                                    u, edges2[u], ugraph.degree(u), mtype2);
                                                // work on update graph only
                                                //if ( reduced_degree2[u] && reduced_degree2[v] )
                                                sorted_intersection(u, red_edges2[u], reduced_degree2[u],
                                                                    v, red_edges2[v], reduced_degree2[v], mtype3);
                                        });
                
                //countUpdateTriangles(ugraph,insertion);
                computeTriangleCount();
        }

        


        // void DynTriangleCounting::updateBatch(const std::vector<GraphEvent>& batch) {

        //         o_t = t_t;
        //         bool insertion = true;
        //         if (batch[0].type == GraphEvent::EDGE_REMOVAL)
        //                 insertion = false;
        //         else if (batch[0].type != GraphEvent::EDGE_ADDITION)
        //                 throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
                
                
        //         Graph ugraph = Graph(G->upperNodeIdBound());
        //         for(auto e : batch){
        //                 // check that all events of the batch are of the same type.
        //                 if (e.type == GraphEvent::EDGE_ADDITION && !insertion) {
        //                         throw std::runtime_error("All Events of an insert batch should be of type EDGE_ADDITION.");
        //                 } else if (e.type == GraphEvent::EDGE_REMOVAL && insertion) {
        //                         throw std::runtime_error("All Events of a remove batch should be of type EDGE_REMOVAL.");

        //                 }
        //                 if (e.type != GraphEvent::EDGE_ADDITION && e.type != GraphEvent::EDGE_REMOVAL) {
        //                         throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
        //                 }


        //                 ugraph.addEdge(e.u, e.v);
        //         }

        //         ugraph.sortEdges();
        //         INFO("** ***** SORTING BATCH *****");

        //         Aux::Timer timer;
        //         if (graphAlreadyUpdated) INFO("** ***** GRAPH ALREADY UPDATED *****");
        //         else {

        //                 timer.start();
        //                 if (insertion) {
        //                         Graph GI = edgeInsertionSorted(ugraph);
        //                         U = std::move(GI);
        //                         INFO(" ** ***** INSERTION RUN SUCCESSFULLY ***** ");
        //                         //INFO(" ** ***** GI  = (" , GI.numberOfNodes() , ", " , GI.numberOfEdges() , ")");
                                
        //                 } else {
        //                         Graph GD = edgeDeletionSorted(ugraph);
        //                         U = std::move(GD);
        //                         INFO(" ** ***** DELETION RUN SUCCESSFULLY ***** ");
        //                         //INFO(" ** ***** GD  = (" , GD.numberOfNodes() , ", " , GD.numberOfEdges() , ")");
        //                 }
        //                 G = &U;
        //                 timer.stop();
        //                 INFO(" ** ***** Time to insert batch =  ", timer.elapsedMicroseconds() / 1e6, " secs.");  
        //         }
        //         // INFO(" ** ***** AFTER GRAPH ASSIGNMENT ***** ");
        //         // INFO(" ** ***** G = (" , G->numberOfNodes() , ", " , G->numberOfEdges() , ")");
        //         // G->forNodes([&](node u) {
        //         //                     G->forEdgesOf(u, [&](node, node v, edgeid) {
        //         //                                             std::cout << " " << v;
        //         //                                     });
        //         //                     std::cout << std::endl;
        //         //             });
                
                
                
        //         countUpdateTriangles(ugraph,insertion);
        //         computeTriangleCount();
        // }


        count DynTriangleCounting::sorted_intersection(node u, const std::vector<node> & u_adj, count u_deg,
                                                       node v, const std::vector<node> & v_adj, count v_deg, double m) {
                count triangles = 0;
                index i = 0, j = 0;
                assert( u_adj.size() == u_deg );
                assert( v_adj.size() == v_deg );

                //assert(u_deg <= u_adj.size());
                //assert(v_deg <= v_adj.size());

                // std::cout << " e: ( " << u << ", " << v << " ) --> list ( [";
                // for (int k = 0; k < u_deg; k++)
                //         std::cout << u_adj[k] << " ";
                // std::cout << "] , [";
                // for (int k = 0; k < v_deg; k++)
                //         std::cout << v_adj[k] << " ";
                // std::cout << "]) "<<  std::endl;


                
                while( (i < u_deg) && (j < v_deg) ) {
                        int comp;
                        comp = u_adj[i] - v_adj[j];
                        triangles += (comp == 0);
                        if (comp == 0) {
#pragma omp atomic
                                TrianglesPerNode[u_adj[i]] += m;
                        }
                        i += (comp <= 0);
                        j += (comp >= 0);
                }
#pragma omp atomic
                TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                TrianglesPerNode[v] += m * triangles;
                return triangles;
        }








        void DynTriangleCounting::countUpdateTriangles(const Graph &ugraph, bool insertion) {
                double mtype1, mtype2, mtype3;
                INFO(" ** ***** UPDATE ");
                INFO(" ** ***** G = (" , G->numberOfNodes() , ", " , G->numberOfEdges() , ")");
                if (insertion) {
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        if (!t_t) return;
                }
                // create local full and reduced copies for G
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                std::vector<std::vector<node> > red_edges(G->upperNodeIdBound());
                std::vector<node> reduced_degree(G->upperNodeIdBound(), 0);
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            red_edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     if (v < u) {
                                                                             red_edges[u].emplace_back(v);
                                                                             reduced_degree[u]++;        
                                                                     }
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });
                // create local full and reduced copies for update graph
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                std::vector<std::vector<node> > red_edges2(ugraph.upperNodeIdBound());
                std::vector<node> reduced_degree2(ugraph.upperNodeIdBound(), 0);
                ugraph.parallelForNodes([&](node u) {
                                                edges2[u].reserve(ugraph.degree(u));
                                                red_edges2[u].reserve(ugraph.degree(u));
                                                ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                             if (v < u) {
                                                                                     red_edges2[u].emplace_back(v);
                                                                                     reduced_degree2[u]++;
                                                                             }
                                                                             edges2[u].emplace_back(v);
                                                                     });
                                        });

                ugraph.parallelForEdges([&](node u, node v) {
                                                
                                                // work on updated graph only
                                                sorted_intersection(u,edges[u], G->degree(u), v, edges[v], G->degree(v), mtype1); 
                                                // work on updated and update graph
                                                sorted_intersection(u, red_edges[u], reduced_degree[u],
                                                                    v, edges2[v], ugraph.degree(v), mtype2);
                                                sorted_intersection(v, red_edges[v], reduced_degree[v],
                                                                    u, edges2[u], ugraph.degree(u), mtype2);
                                                // work on update graph only
                                                //if ( reduced_degree2[u] && reduced_degree2[v] )
                                                sorted_intersection(u, red_edges2[u], reduced_degree2[u],
                                                                    v, red_edges2[v], reduced_degree2[v], mtype3);
                                        });


        }



        // deletion merge for a node based on adjacency list in G and B
        count DynTriangleCounting::do_merge_delete( std::vector<node> & adj_inG, const count deg_inG,
                                                   const std::vector<node> & adj_inB, const count deg_inB) {
                assert(adj_inG.size() >= deg_inG);
                assert(adj_inB.size() >= deg_inB);
                
                index a = 0, b = 0, k = 0;
                while (a < deg_inG && b < deg_inB) {
                        int comp = adj_inG[a] - adj_inB[b];
                        if (!comp) { // deleted
                                a++; b++;
                        } else if ( comp < 0 ) {
                                adj_inG[k] = adj_inG[a];
                                k++; a++;
                        } else { // deletion not found
                                b++;
                        }
                }
                while (a < deg_inG) {
                        adj_inG[k] = adj_inG[a];
                        k++; a++;
                }
                while (b < deg_inB) {
                        b++;
                }
                index r = k;
                // erasing unnecessary elements
                while (r < deg_inG) {
                        adj_inG.pop_back();
                        r++;
                }
                
                return k;
        }
        

        
        // insertion merge for a node based on adjacency list in G and B
        count DynTriangleCounting::do_merge_insert(const std::vector<node> adj_inG, const count deg_inG,
                                                   const std::vector<node> adj_inB, const count deg_inB,
                                                   std::vector<node> & adj_inGI) {
                assert(adj_inG.size() >= deg_inG);
                assert(adj_inB.size() >= deg_inB);
                assert(adj_inGI.size() <= (deg_inG + deg_inB));
                assert(adj_inGI.size() >= deg_inG);
                
                index a = 0, b = 0, k = 0;
                while (a < deg_inG && b < deg_inB) {
                        int comp = adj_inG[a] - adj_inB[b];
                        if (!comp) { // this edge is already present in graph.
                                adj_inGI[k] = adj_inG[a];
                                k++; a++; b++;
                        } else if ( comp < 0 ) {
                                adj_inGI[k] = adj_inG[a];
                                k++; a++;
                        } else {
                                adj_inGI[k] = adj_inB[b];
                                k++; b++;
                        }
                }
                while (a < deg_inG) {
                        adj_inGI[k] = adj_inG[a];
                        k++; a++;
                }
                while (b < deg_inB) {
                        adj_inGI[k] = adj_inB[b];
                        k++; b++;
                }
                return k;
        }
        
        Graph  DynTriangleCounting::edgeInsertionSorted(const Graph & ugraph) {
                count n = G->upperNodeIdBound();
                assert(checkSorted(NULL));
                assert(checkSorted(&ugraph));

                std::vector<edgeweight> empty;
                GraphBuilder result(n, false, false);
                
                //INFO(" ** ***** G = (" , G->numberOfNodes() , ", " , G->numberOfEdges() , ")");
                G->balancedParallelForNodes([&](node v) {

                                                    count current_degree = G->degree(v);
                                                    count update_degree = ugraph.degree(v);
                                                    count max_new_degree = current_degree + update_degree;
                                                    // allocate
                                                    std::vector<node> ngb;
                                                    std::vector<node> update_ngb;
                                                    ngb.reserve(current_degree);
                                                             
                                                    G->forNeighborsOf(v, [&](node _i, node j) {
                                                                                 ngb.emplace_back(j);
                                                                         });
                                                    
                                                    // INFO("** ***** NODE v = ", v);
                                                    // INFO("** ***** NODE degree = ", current_degree);
                                                    // std::cout << "[";
                                                    // for (int i = 0; i < ngb.size(); i++) {
                                                    //         std::cout << ngb[i] << " ";
                                                    // }
                                                    // std::cout << "]" << std::endl;

                                                    // TODO: should check also if the other degree is zero to avoid comparisons !
                                                    
                                                    if (update_degree) {
                                                            update_ngb.reserve(update_degree);
                                                            
                                                            

                                                            ugraph.forNeighborsOf(v, [&](node _i, node j) {
                                                                                             update_ngb.emplace_back(j);
                                                                                     });
                                                            
                                                            std::vector<node> new_ngb(max_new_degree); // maybe reserve ?

                                                            // INFO("** ***** MERGING ", v);
                                                            // std::cout << "[";
                                                            // for (int i = 0; i < update_ngb.size(); i++) {
                                                            //         std::cout << update_ngb[i] << " ";
                                                            // }
                                                            // std::cout << "]" << std::endl;
                                                            
                                                            count new_degree = do_merge_insert(ngb, current_degree, update_ngb, update_degree, new_ngb);
                                                            // INFO("** ***** AFTER MERGING OF ", v);
                                                            // std::cout << "[";
                                                            // for (int i = 0; i < new_ngb.size(); i++) {
                                                            //         std::cout << new_ngb[i] << " ";
                                                            // }
                                                            // std::cout << "]" << std::endl;
                                                            
                                                            assert(new_ngb.size() == new_degree);


                                                            result.swapNeighborhood(v, new_ngb, empty, false);
                                                    } // just copy previous adj
                                                    
                                                    else {
                                                            // INFO("** ***** NO MERGING ", v);
                                                            // std::cout << "[";
                                                            // for (int i = 0; i < ngb.size(); i++) {
                                                            //         std::cout << ngb[i] << " ";
                                                            // }
                                                            // std::cout << "]" << std::endl;
                                                            result.swapNeighborhood(v, ngb, empty, false);
                                                    }
                                            });
                return result.toGraph(false,true);
                                                            
        }




                Graph  DynTriangleCounting::edgeDeletionSorted(const Graph & ugraph) {
                count n = G->upperNodeIdBound();
                assert(checkSorted(NULL));
                assert(checkSorted(&ugraph));

                std::vector<edgeweight> empty;
                GraphBuilder result(n, false, false);

                //INFO(" ** ***** G = (" , G->numberOfNodes() , ", " , G->numberOfEdges() , ")");
                G->balancedParallelForNodes([&](node v) {

                                                    count current_degree = G->degree(v);
                                                    count update_degree = ugraph.degree(v);
                                                    count max_new_degree = current_degree + update_degree;
                                                    // allocate
                                                    std::vector<node> ngb;
                                                    std::vector<node> update_ngb;
                                                    ngb.reserve(current_degree);
                                                             
                                                    G->forNeighborsOf(v, [&](node _i, node j) {
                                                                                 ngb.emplace_back(j);
                                                                         });
                                                    
                                                    // INFO("** ***** NODE v = ", v);
                                                    // INFO("** ***** NODE degree = ", current_degree);
                                                    // std::cout << "[";
                                                    // for (int i = 0; i < ngb.size(); i++) {
                                                    //         std::cout << ngb[i] << " ";
                                                    // }
                                                    // std::cout << "]" << std::endl;

                                                   
                                                    if (update_degree) {
                                                            update_ngb.reserve(update_degree);
                                                            ugraph.forNeighborsOf(v, [&](node _i, node j) {
                                                                                             update_ngb.emplace_back(j);
                                                                                     });
                                                            
                                                            // INFO("** ***** MERGING ", v);
                                                            // std::cout << "[";
                                                            // for (int i = 0; i < update_ngb.size(); i++) {
                                                            //         std::cout << update_ngb[i] << " ";
                                                            // }
                                                            // std::cout << "]" << std::endl;
                                                            
                                                            count new_degree = do_merge_delete(ngb, current_degree, update_ngb, update_degree);
                                                            // INFO("** ***** AFTER MERGING OF ", v);
                                                            // std::cout << "[";
                                                            // for (int i = 0; i < ngb.size(); i++) {
                                                            //         std::cout << ngb[i] << " ";
                                                            // }
                                                            // std::cout << "]" << std::endl;
                                                            assert(new_degree == ngb.size());
                                                    }
                                                    // just copy previous adj
                                                    // better move instead of swap
                                                    result.swapNeighborhood(v, ngb, empty, false);

                                            });
                return result.toGraph(false,true);
                                                            
        }
        
        void DynTriangleCounting::edgeInsertion(const std::vector<GraphEvent>& batch) {
                for(auto e : batch){
                        G->addEdge(e.u, e.v);
                }
                if(!checkSorted(NULL)) {
                        G->sortEdges();
                        INFO("** ***** SORTING GRAPH *****");
                }
                else {
                        INFO("** ***** NO NEED TO SORT AFTER INSERTION *****");
                }
        }

        void DynTriangleCounting::edgeDeletion(const std::vector<GraphEvent>& batch) {
                for(auto e : batch){
                        G->removeEdge(e.u, e.v);
                }
                if(!checkSorted(NULL)) {
                        G->sortEdges();
                        INFO("** ***** SORTING GRAPH *****");
                }
                else {
                        INFO("** ***** NO NEED TO SORT AFTER DELETION *****");
                }
        }

        
        bool DynTriangleCounting::checkSorted(const Graph * G1 = NULL) {
                bool res = true;
                if (G1 == NULL) {
                        // not breaking involved
                        G->forNodes([&] (node u) {
                                           if(!std::is_sorted(G->neighborRange(u).begin(), G->neighborRange(u).end()))
                                                   res = false;
                                    });
                }
                else {
                        G1->forNodes([&] (node u) {
                                            if(!std::is_sorted(G1->neighborRange(u).begin(), G1->neighborRange(u).end()))
                                                    res = false;
                                     });
                }
                return res;
        }
                                   
                        



        
} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_
