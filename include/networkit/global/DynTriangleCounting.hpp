/*
 * DynKatzCentrality.hpp
 *
 *  Created on: April 2018
 *      Author: Alexander van der Grinten
 *      based on code by Elisabetta Bergamini
 */

#ifndef NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_
#define NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_

#include <omp.h>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>

namespace NetworKit {

/**
 * Calculates the number of triangles in a dynamic graph.
 *
 * @ingroup global
 */
class DynTriangleCounting : public Algorithm, public DynAlgorithm {

public:
  /**
   * Constructs the DynTopHarmonicCloseness class. This class implements dynamic
   * algorithms for harmonic top-k closeness centrality. The implementation
   * is based on the static algorithms by Borassi et al. (complex networks)
   * and Bergamini et al. (large-diameter networks).
   *
   * @param G The graph.
   * @param k The number of top nodes.
   * @param useBFSbound Whether to use the algorithm for networks with large
   * diameter
   */
        
        DynTriangleCounting(Graph &G, bool insertion = true)
                : G(&G), insertion(insertion),
                  TrianglesPerNode(G.upperNodeIdBound(), 0.0),
                  t_t(0.0), S1(0), S2(0), S3(0) {
                std::cout << " TC: set up for insertion (by default). " << std::endl;
        }
        
        void run() override;
        
        /* insertion function works as a reset for insert/deletion. */
        void reset(Graph &G, bool insert = true) {
                insertion = insert;
                if (insertion) std::cout << " TC : reset for insertion. " << std::endl;
                else std::cout << " TC : reset for deletion. " << std::endl;
                std::fill(TrianglesPerNode.begin(), TrianglesPerNode.end(), 0.0);
                TrianglesPerNode.resize(G.upperNodeIdBound(), 0.0);
                // TODO: check, clear only for insertion. Not for removal
                //if (insertion) t_t = 0;
                t_t = 0.0;
                S1 = S2 = S3 = 0; 
        }
        // tmp version for inserting edges and sorting afterwards
        void edgeInsertion(const std::vector<GraphEvent> &batch);
        void edgeDeletion(const std::vector<GraphEvent> &batch);
        void edgeInsertionSorted(const std::vector<GraphEvent> &batch);
        void edgeDeletionSorted(const std::vector<GraphEvent> &batch);
        
        
        bool Insertion() {
                return insertion;
        }
        /**
         * Updates the list of k nodes with the highest closeness in G
         * after a batch of updates.
         *
         * @param batch A vector of edge modification events.
         */
        void updateBatch(const std::vector<GraphEvent> &batch) override;


        /**
         * Updates the list of the k nodes with the highest closeness in G.
         *
         * @param event The edge modification event.
         */
        
        void update(GraphEvent singleEvent) override {
                std::vector<GraphEvent> events{singleEvent};
                updateBatch(events);
        }

        /**
         * Get number of triangles
         */
        count getTriangleCount() {
                assureFinished();
                return t_t;
        }
        count computeTriangleCount() {
                t_t = G->parallelSumForNodes([&](node u) {
                                                     return TrianglesPerNode[u];
                                             });
                
                t_t = t_t/3;
        }
        

        std::vector<double> getTriangleScores() {
                assureFinished();
                return TrianglesPerNode;
        }
        bool checkSorted(const Graph * G1);

        
private:

        count bsearch_intersection( node u, node v, const std::vector<node> & u_adj, count u_deg);
        count sorted_intersection( node u, const std::vector<node> & u_adj, count u_deg,
                                   node v,const std::vector<node> & v_adj, count v_deg, double m);
        
        count countTriangleN1O2(const Graph & G, double m);
        count countTriangleN2O1(const Graph & G, double m);
        count countTriangleN3(const Graph & G, double m);

        count countTrianglesType1(const Graph &ugraph, double m);
        count countTrianglesType2(const Graph &ugraph, double m);
        count countTrianglesType3(const Graph &ugraph, double m);
        void countUpdateTriangles(const Graph &ugraph);


        Graph *G;
        bool insertion;
        std::vector<double> TrianglesPerNode;
        double t_t;
        count S1, S2, S3;
        

        

};

        // form stringer 
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
                
                std::cout << "TriCnt : RUN " << std::endl;
                std::fill(TrianglesPerNode.begin(),  TrianglesPerNode.end(), 0.0);
                // resize in case of change in graph size
                TrianglesPerNode.resize(G->upperNodeIdBound(), 0.0);
        
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                /* sort edge lists */
                // TODO: remove it from here. Assume that graph input is sorted
                G->sortEdges();
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
        }

        // update routine.
        void DynTriangleCounting::updateBatch(const std::vector<GraphEvent>& batch) {
                
                count n = G->upperNodeIdBound();
                Graph ugraph = Graph(n);
                for(auto e : batch){
                        ugraph.addEdge(e.u, e.v);
                }
                ugraph.sortEdges();
                assert(checkSorted(&ugraph));

                double mtype1, mtype2, mtype3;
                int total = 0;
                int S1 = 0;
                int S2 = 0;
                int S3 = 0;
                // setting multiplicative parameters for each type of triangle and each mode: insertion/deletion
                if (insertion) {
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                        std::cout << "TriCnt : UPDATE TO INSERT - PREVIOUS COUNT =  " << t_t << std::endl;

                        countTrianglesType1(ugraph, mtype1);
                        countTrianglesType2(ugraph, mtype2);
                        countTrianglesType3(ugraph, mtype3);
                        //t_t += (S1-S2+S3);
                        
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        std::cout << "TriCnt : UPDATE TO DELETE - PREVIOUS COUNT =  " << t_t << std::endl;
                        if (!t_t) return;
                        countTrianglesType1(ugraph, mtype1);
                        countTrianglesType2(ugraph, mtype2);
                        countTrianglesType3(ugraph, mtype3);
                        // t_t = S1+S2+S3
                }


                //countUpdateTriangles(ugraph);
                computeTriangleCount();
                std::cout << "TriCnt : FINAL COUNT =  " << t_t << std::endl;
        }


        count DynTriangleCounting::sorted_intersection(node u, const std::vector<node> & u_adj, count u_deg,
                                                       node v, const std::vector<node> & v_adj, count v_deg, double m) {
                count triangles = 0;
                index i = 0, j = 0;
                assert(u_adj.size() <= u_deg);
                assert(v_adj.size() <= v_deg);
                
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








        void DynTriangleCounting::countUpdateTriangles(const Graph &ugraph) {
                double mtype1, mtype2, mtype3;
                int total = 0;
                int S1 = 0;
                int S2 = 0;
                int S3 = 0;
                // setting multiplicative parameters for each type of triangle and each mode: insertion/deletion
                if (insertion) {
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                        std::cout << "TriCnt : UPDATE TO INSERT - PREVIOUS COUNT =  " << t_t << std::endl;
                        //t_t += (S1-S2+S3);
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        std::cout << "TriCnt : UPDATE TO DELETE - PREVIOUS COUNT =  " << t_t << std::endl;
                        if (!t_t) return;
                        // t_t = S1+S2+S3
                }
                // create local copies and reduced degrees for G
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                std::vector<node> reduced_degree(G->upperNodeIdBound(), 0);
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     if (v < u) {
                                                                             edges[u].emplace_back(v);
                                                                             reduced_degree[u]++;        
                                                                     }
                                                             });
                                    });
                // create local copies and reduced degrees for update graph ugraph
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                std::vector<node> reduced_degree2(ugraph.upperNodeIdBound(), 0);
                ugraph.parallelForNodes([&](node u) {
                                                edges2[u].reserve(ugraph.degree(u));
                                                ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                             if (v < u) {
                                                                                     edges2[u].emplace_back(v);
                                                                                     reduced_degree2[u]++;
                                                                             }
                                                                     });
                                        });
                
                ugraph.parallelForEdges([&](node u, node v) {
                                                // work on updated graph only
                                                S1 = sorted_intersection(u,edges[u], G->degree(u), v, edges[v], G->degree(v), mtype1); 
                                                // work on updated and update graph
                                                S2 = sorted_intersection(u, edges[u], reduced_degree[u],
                                                                    v, edges2[v], ugraph.degree(v), mtype2); 
                                                S2 = S2 + sorted_intersection(v, edges[v], reduced_degree[v],
                                                                    u, edges2[u], ugraph.degree(u), mtype2);
                                                // work on update graph only
                                                S3 = sorted_intersection(u, edges2[u], reduced_degree2[u],
                                                                    v, edges2[v], reduced_degree2[v], mtype3);
                                        });    

                if (insertion) {
                        t_t += S1-S2+S3;
                        std::cout << "TriCnt : FINAL =  " << t_t << std::endl;
                } else {
                        t_t = S1+S2+S3;
                        std::cout << "TriCnt : FINAL =  " << t_t << std::endl;
                } 
        }

        
        

        count DynTriangleCounting::countTrianglesType1(const Graph &ugraph, double m) {
                count total =  0;
                // following is repetiton compared to countTriangleN2O1 and countTriangleN3
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });
                
                ugraph.parallelForEdges([&](node u, node v) {
                                           total += sorted_intersection(u,edges[u], G->degree(u), v, edges[v], G->degree(v), m);
                                           
                                   });
                return total;
        }

        count DynTriangleCounting::countTrianglesType2(const Graph &ugraph, double m) {
                count total =  0;
                std::vector<node> reduced_degree(G->upperNodeIdBound(), 0);
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     if (v < u) {
                                                                             edges[u].emplace_back(v);
                                                                             reduced_degree[u]++;        
                                                                     }
                                                             });
                                    });
                
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                ugraph.parallelForNodes([&](node u) {
                                           edges2[u].reserve(ugraph.degree(u));
                                           ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                   edges2[u].emplace_back(v);
                                                           });
                                   });

                ugraph.parallelForEdges([&](node u, node v) {                           
                                                  count d_u = reduced_degree[u];
                                                  total += sorted_intersection(u, edges[u], d_u, v, edges2[v], ugraph.degree(v), m); 
                                                  count d_v = reduced_degree[v];
                                                  total += sorted_intersection(v, edges[v], d_v, u, edges2[u], ugraph.degree(u), m); 

                                          });
                return total;
        }

                count DynTriangleCounting::countTrianglesType3(const Graph &ugraph, double m) {
                count total =  0;
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                std::vector<node> reduced_degree(ugraph.upperNodeIdBound(), 0);
                
                ugraph.parallelForNodes([&](node u) {
                                                  edges2[u].reserve(ugraph.degree(u));
                                                  ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                                 if (v < u) {
                                                                                         edges2[u].emplace_back(v);
                                                                                         reduced_degree[u]++;
                                                                                 }
                                                                         });
                                          });
                
                ugraph.parallelForEdges([&](node u, node v) {
                                           
                                           count d_u = reduced_degree[u];
                                           count d_v = reduced_degree[v];
                                           total += sorted_intersection(u, edges2[u], d_u, v, edges2[v], d_v, m);
                                   });
                return total;
        }



        

        

//         count DynTriangleCounting::countTriangleN1O2(const Graph &B, double m) {
//                 count total =  0;
//                 // following is repetiton compared to countTriangleN2O1 and countTriangleN3
//                 std::vector<std::vector<node> > edges(G->upperNodeIdBound());
//                 G->parallelForNodes([&](node u) {
//                                             edges[u].reserve(G->degree(u));
//                                             G->forEdgesOf(u, [&](node, node v, edgeid) {
//                                                                      edges[u].emplace_back(v);
//                                                              });
//                                     });
                
//                 B.parallelForEdges([&](node u, node v) {
//                                            count triangles = 0;
//                                            count d_u = G->degree(u); 
//                                            count d_v = G->degree(v); 
//                                            index i = 0, j = 0;
//                                            std::cout << "TC N1O2 : Edge ( " << u << ", " << v  << " )"  << std::endl;
//                                            while( (i < d_u) && (j < d_v) ) {
//                                                    int comp;
//                                                    comp = edges[u][i] - edges[v][j];
//                                                    triangles += (comp == 0);
//                                                    if (comp == 0) {
// #pragma omp atomic
//                                                            TrianglesPerNode[edges[u][i]] += m;
//                                                            std::cout << "TC N1O2 : triangle count = " << triangles << "( due to common neighbor " << edges[u][i] << ")" << std::endl;
//                                                    }
//                                                    i += (comp <= 0);
//                                                    j += (comp >= 0);
//                                            }
// #pragma omp atomic
//                                            TrianglesPerNode[u] += m * triangles;
// #pragma omp atomic
//                                            TrianglesPerNode[v] += m * triangles;

//                                            total += 3* triangles; 
//                                    });
//                 std::cout << "TC N1O2 : total = " << total << std::endl;
                
//                 return total;
//         }



        count DynTriangleCounting::countTriangleN1O2(const Graph &B, double m) {
                count total =  0;
                std::cout << "TriCnt : TYPE DELTA1 TRIANGLES " << std::endl;
                std::cout << "TriCnt : PARAMETER M = " << m  << std::endl;
                
                // following is repetiton compared to countTriangleN2O1 and countTriangleN3
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });

                
                B.parallelForEdges([&](node u, node v) {
                                           count triangles = 0;
                                           count d_u = G->degree(u); 
                                           count d_v = G->degree(v); 
                                           index i = 0, j = 0;
                                           //std::cout << "TC N1O2 : Edge ( " << u << ", " << v  << " )"  << std::endl;
                                           while( (i < d_u) && (j < d_v) ) {
                                                   int comp;
                                                   comp = edges[u][i] - edges[v][j];
                                                   triangles += (comp == 0);
                                                   if (comp == 0) {
#pragma omp atomic
                                                           TrianglesPerNode[edges[u][i]] +=  m;
                                                           //std::cout << "TC N1O2 : triangle count = " << triangles << "( due to common neighbor w = " << edges[u][i] << ") -- TrianglesPerNode[w] = " << TrianglesPerNode[edges[u][i]]  << " m = " << m << std::endl;
                                                           
                                                   }
                                                   i += (comp <= 0);
                                                   j += (comp >= 0);
                                           }
#pragma omp atomic
                                           TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                                           TrianglesPerNode[v] += m * triangles;

                                           total += triangles; 
                                   });
                return total;
        }



        count DynTriangleCounting::countTriangleN2O1(const Graph &B, double m) {
                count total =  0;

                std::cout << "TriCnt : TYPE DELTA2 TRIANGLES " << std::endl;
                std::cout << "TriCnt : PARAMETER M = " << m  << std::endl;
                
                
                std::vector<node> reduced_degree(G->upperNodeIdBound(), 0);
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     if (v < u) {
                                                                             edges[u].emplace_back(v);
                                                                             reduced_degree[u]++;        
                                                                     }
                                                             });
                                    });
                
                std::vector<std::vector<node> > edges2(B.upperNodeIdBound());
                B.parallelForNodes([&](node u) {
                                           edges2[u].reserve(B.degree(u));
                                           B.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                   edges2[u].emplace_back(v);
                                                           });
                                   });
                
                B.parallelForEdges([&](node x, node y) {

                                           
                                           node u = x;
                                           node v = y;
                                           // std::cout << "TC N2O1 : Edge ( " << u << ", " << v  << " )" << " degrees ( "
                                           //           << G->degree(u) << ", " << B.degree(v)  << " )"
                                           //           << " reduced_degree for G-HAT ( " << reduced_degree[u] << "  )"<< std::endl;
                                           count triangles = 0;

                                           count d_u = reduced_degree[u];
                                           count d_v = B.degree(v); 
                                           index i = 0, j = 0;
                                           
                                           while( (i < d_u) && (j < d_v) ) {
                                                   int comp;                                                   
                                                   comp = edges[u][i] - edges2[v][j];
                                                   triangles += (comp == 0);
                                                   if (comp == 0) {
#pragma omp atomic
                                                           TrianglesPerNode[edges[u][i]] += m;
                                                           //std::cout << "TC N1O2 : triangle count = " << triangles << "( due to common neighbor " << edges[u][i] << ")" << std::endl;
                                                   }
                                                           // atomicSub(TrianglesPerNode + edges[u][i], 1);
                                                   i += (comp <= 0);
                                                   j += (comp >= 0);
                                           }
#pragma omp atomic
                                           TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                                           TrianglesPerNode[v] += m * triangles;
                                           //atomicAdd(TrianglesPerNode + u, m * triangles);
                                           //atomicAdd(TrianglesPerNode + v, m * triangles);
                                           total +=   triangles; 
                                           //std::cout << "TC N2O1 : (one direction -->) total = " << total << std::endl;

                                           u = y;
                                           v = x;
                                           // std::cout << "TC N2O1 : Edge ( " << u << ", " << v  << " )" << " degrees ( "
                                           //           << G->degree(u) << ", " << B.degree(v)  << " )"
                                           //           << " reduced_degree for G-HAT ( " << reduced_degree[u] << "  )"<< std::endl;
                                           triangles = 0;
                                           
                                           d_u = reduced_degree[u];
                                           d_v = B.degree(v); 
                                           i = 0, j = 0;
                                           
                                           while( (i < d_u) && (j < d_v) ) {
                                                   int comp;                                                   
                                                   comp = edges[u][i] - edges2[v][j];
                                                   triangles += (comp == 0);
                                                   if (comp == 0) {
#pragma omp atomic
                                                           TrianglesPerNode[edges[u][i]] += m;
                                                           //std::cout << "TC N1O2 : triangle count = " << triangles << "( due to common neighbor " << edges[u][i] << ")" << std::endl;
                                                   }
                                                           // atomicSub(TrianglesPerNode + edges[u][i], 1);
                                                   i += (comp <= 0);
                                                   j += (comp >= 0);
                                           }
#pragma omp atomic
                                           TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                                           TrianglesPerNode[v] += m * triangles;
                                           //atomicAdd(TrianglesPerNode + u, m * triangles);
                                           //atomicAdd(TrianglesPerNode + v, m * triangles);
                                           total +=  triangles; 


                                           //std::cout << "TC N2O1 : (other direction <--) total = " << total << std::endl;


                                   });
                return total;
        }


        

        
//         count DynTriangleCounting::countTriangleN2O1(const Graph &B, double m) {
//                 count total =  0;
//                 std::vector<std::vector<node> > edges(G->upperNodeIdBound());
//                 G->parallelForNodes([&](node u) {
//                                             edges[u].reserve(G->degree(u));
//                                             G->forEdgesOf(u, [&](node, node v, edgeid) {
//                                                                      edges[u].emplace_back(v);
//                                                              });
//                                     });
                
//                 std::vector<std::vector<node> > edges2(B.upperNodeIdBound());
//                 B.parallelForNodes([&](node u) {
//                                            edges2[u].reserve(B.degree(u));
//                                            B.forEdgesOf(u, [&](node, node v, edgeid) {
//                                                                    edges2[u].emplace_back(v);
//                                                            });
//                                    });
                
//                 B.parallelForEdges([&](node u, node v) {
//                                            count triangles = 0;
//                                            count d_u = G->degree(u); // G or B ?
//                                            count d_v = B.degree(v); // G or B ?
//                                            index i = 0, j = 0;
                                           
//                                            while( (i < d_u) && (j < d_v) ) {
//                                                    int comp;                                                   
//                                                    comp = edges[u][i] - edges2[v][j];
//                                                    triangles += (comp == 0);
//                                                    if (comp == 0)
// #pragma omp atomic
//                                                            TrianglesPerNode[edges[u][i]] += m;
//                                                            // atomicSub(TrianglesPerNode + edges[u][i], 1);
//                                                    i += (comp <= 0);
//                                                    j += (comp >= 0);
//                                            }
// #pragma omp atomic
//                                            TrianglesPerNode[u] += m * triangles;
// #pragma omp atomic
//                                            TrianglesPerNode[v] += m * triangles;
//                                            //atomicAdd(TrianglesPerNode + u, m * triangles);
//                                            //atomicAdd(TrianglesPerNode + v, m * triangles);
//                                            total += 3* triangles; 
//                                    });
//                 return total;
//         }



                count DynTriangleCounting::countTriangleN3(const Graph &B, double m) {

                        std::cout << "TriCnt : TYPE DELTA3 TRIANGLES " << std::endl;
                        std::cout << "TriCnt : PARAMETER M = " << m  << std::endl;
                        count total =  0;
                
                        std::vector<std::vector<node> > edges2(B.upperNodeIdBound());
                        std::vector<node> reduced_degree(B.upperNodeIdBound(), 0);
                        std::cout << "TC N3 : m = " << m  << std::endl;
                        B.parallelForNodes([&](node u) {
                                           edges2[u].reserve(B.degree(u));
                                           B.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                   if (v < u) {
                                                                   edges2[u].emplace_back(v);
                                                                   reduced_degree[u]++;
                                                                   }
                                                           });
                                   });
                
                B.parallelForEdges([&](node u, node v) {
                                           count triangles = 0;
                                           count d_u = reduced_degree[u];
                                           count d_v = reduced_degree[v];
                                           index i = 0, j = 0;
                                           // std::cout << "TC N3 : Edge ( " << u << ", " << v  << " )" << " degrees ( "
                                           //           << B.degree(u) << ", " << B.degree(v)  << " )"
                                           //           << " reduced_degree for G' ( " << reduced_degree[u] << ", "
                                           //           << reduced_degree[v] << "  )"<< std::endl;

                                           
                                           while( (i < d_u) && (j < d_v) ) {
                                                   int comp;                                                   
                                                   comp = edges2[u][i] - edges2[v][j];
                                                   triangles += (comp == 0);
                                                   if (comp == 0) {
#pragma omp atomic
                                                           TrianglesPerNode[edges2[u][i]] += m;
                                                           //std::cout << "TC N3 : triangle count = " << triangles << "( due to common neighbor " << edges2[u][i] << ")" << std::endl;
                                                   }
                                                   i += (comp <= 0);
                                                   j += (comp >= 0);
                                           }
#pragma omp atomic
                                           TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                                           TrianglesPerNode[v] += m * triangles;
                                           total +=  triangles; 
                                   });
                return total;
        }



        
     
        

        
//         count DynTriangleCounting::countTriangleN3(const Graph &B, double m) {
//                 count total =  0;
//                 std::vector<std::vector<node> > edges2(B.upperNodeIdBound());
//                 B.parallelForNodes([&](node u) {
//                                            edges2[u].reserve(B.degree(u));
//                                            B.forEdgesOf(u, [&](node, node v, edgeid) {
//                                                                    edges2[u].emplace_back(v);
//                                                            });
//                                    });
                
//                 B.parallelForEdges([&](node u, node v) {
//                                            count triangles = 0;
//                                            count d_u = B.degree(u); // G or B ?
//                                            count d_v = B.degree(v); // G or B ?
//                                            index i = 0, j = 0;
                                           
//                                            while( (i < d_u) && (j < d_v) ) {
//                                                    int comp;                                                   
//                                                    comp = edges2[u][i] - edges2[v][j];
//                                                    triangles += (comp == 0);
//                                                    if (comp == 0)
// #pragma omp atomic
//                                                            TrianglesPerNode[edges2[u][i]] += m;
//                                                    i += (comp <= 0);
//                                                    j += (comp >= 0);
//                                            }
// #pragma omp atomic
//                                            TrianglesPerNode[u] += m * triangles;
// #pragma omp atomic
//                                            TrianglesPerNode[v] += m * triangles;
//                                            total += 3* triangles; 
//                                    });
//                 return total;
//         }


        


        

        


        void DynTriangleCounting::edgeInsertionSorted(const std::vector<GraphEvent>& batch) {
                
                // construct G'
                count n = G->upperNodeIdBound();
                Graph B = Graph(n);
                for(auto e : batch){
                        B.addEdge(e.u, e.v);
                }
                B.sortEdges();
                if(checkSorted(NULL)) std::cout << " TC : original graph G is sorted. " << std::endl;
                if(checkSorted(&B)) std::cout << " TC : batch-graph B is sorted. "<< std::endl;
                
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                std::vector<std::vector<node> > edges2(B.upperNodeIdBound());

                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });

                B.parallelForNodes([&](node u) {
                                            edges2[u].reserve(B.degree(u));
                                            B.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                    edges2[u].emplace_back(v);
                                                            });
                                   });

                std::cout << " TC : retrieving edge lists. "<< std::endl;
                G->parallelForNodes([&](node u) {
                                            int i = G->degree(u);
                                            int j = B.degree(u);

                                            while( (i >= 0) && (j >= 0) ) {
                                                    int comp;
                                                    node x = G->upperNodeIdBound();
                                                    node y = G->upperNodeIdBound();
                                                    if ( i == 0 && G->isIsolated(u)) x = 0;  
                                                    else x = edges[u][i];
                                                    if ( j == 0 && B.isIsolated(u)) y = 0;  
                                                    else y = edges2[u][j];
                                                    comp = x - y;
                                                    //comp = edges[u][i] - edges2[u][j];
                                                   if (comp > 0) {

                                                           G->addEdge(u, y);
                                                           i -= 1; 
                                                   } else {
                                                           G->addEdge(u, edges2[u][j]);
                                                           j -= 1; 
                                                   }
                                           }
                                            while(j > 0) {
                                            //while( j >= 0) {
                                                   G->addEdge(u, edges2[u][j]);
                                                   j -= 1; 
                                           }
                                    });
                std::cout << " TC : merging sorted graphs sucessfully. "<< std::endl;
                
        }


        void DynTriangleCounting::edgeInsertion(const std::vector<GraphEvent>& batch) {
                for(auto e : batch){
                        G->addEdge(e.u, e.v);
                }
                G->sortEdges();
                std::cout << " TC : sorting again after each batch insertion. "<< std::endl;
                
        }

        void DynTriangleCounting::edgeDeletion(const std::vector<GraphEvent>& batch) {
                for(auto e : batch){
                        G->removeEdge(e.u, e.v);
                }
                G->sortEdges();
                std::cout << " TC : sorting again after each batch deletion. "<< std::endl;
                
        }


        
        void DynTriangleCounting::edgeDeletionSorted(const std::vector<GraphEvent> &batch) {
                std::cout << " TC: edgeDeletionSorted() is not implemented yet." << std::endl;
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