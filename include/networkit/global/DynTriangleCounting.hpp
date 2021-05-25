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
        
        DynTriangleCounting(Graph &G)
                : G(&G), TrianglesPerNode(G.upperNodeIdBound(), 0.0),
                  t_t(0.0), o_t(0.0) {
        }
        
        void run() override;

        

        // tmp version for inserting edges and sorting afterwards
        void edgeInsertion(const std::vector<GraphEvent> &batch);
        void edgeDeletion(const std::vector<GraphEvent> &batch);
        void edgeInsertionSorted(const std::vector<GraphEvent> &batch);
        void edgeDeletionSorted(const std::vector<GraphEvent> &batch);
        
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
         * Get number of triangles
         */
        count getTriangleCount() {
                assureFinished();
                return t_t;
        }

        count getNewTriangles() {
                assureFinished();
                return abs(t_t-o_t);
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
        

        count countTrianglesType1(const Graph &ugraph, double m);
        count countTrianglesType2(const Graph &ugraph, double m);
        count countTrianglesType3(const Graph &ugraph, double m);
        void countUpdateTriangles(const Graph &ugraph, bool insertion);


        Graph *G;
        // current triangle count per node
        std::vector<double> TrianglesPerNode;
        // current total number of triangles
        double t_t;
        // old total number of triangles
        double o_t;

        

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
                o_t = t_t;
        }


        void DynTriangleCounting::updateBatch(const std::vector<GraphEvent>& batch) {
                // setting old count equal to current
                o_t = t_t;
                bool insertion = true;
                if (batch[0].type == GraphEvent::EDGE_REMOVAL)
                        insertion = false;
                else if (batch[0].type != GraphEvent::EDGE_ADDITION)
                        throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
                
                
                // create update graph
                Graph ugraph = Graph(G->upperNodeIdBound()); // update graph same size as original one
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

                double mtype1, mtype2, mtype3;
                // setting multiplicative parameters for each type of triangle and each mode: insertion/deletion
                if (insertion) {
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                        countTrianglesType1(ugraph, mtype1);
                        countTrianglesType2(ugraph, mtype2);
                        countTrianglesType3(ugraph, mtype3);
                        //t_t += (S1-S2+S3);
                        
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        if (!t_t) return;
                        count S1 =  countTrianglesType1(ugraph, mtype1);
                        count S2 = countTrianglesType2(ugraph, mtype2);
                        count S3 = countTrianglesType3(ugraph, mtype3);
                        //t_t = S1+S2+S3
                        //std::cout << " total =  " << S1+S2+S3  << std::endl;
                }
                
                //countUpdateTriangles(ugraph,true);
                computeTriangleCount();
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








        void DynTriangleCounting::countUpdateTriangles(const Graph &ugraph, bool insertion) {
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
                }
                else {
                        mtype1 = mtype2 = mtype3 = -1.0;
                        std::cout << "TriCnt : UPDATE TO DELETE - PREVIOUS COUNT =  " << t_t << std::endl;
                        if (!t_t) return;

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
                // following is repetiton compared to countTriangleType2 and countTriangleType3
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->parallelForNodes([&](node u) {
                                            edges[u].reserve(G->degree(u));
                                            G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                     edges[u].emplace_back(v);
                                                             });
                                    });
                
                ugraph.parallelForEdges([&](node u, node v) {
                                                //std::cout << " T1 e: ( " << u << ", " << v << " )" << std::endl;
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
                                                //std::cout << " T2 e: ( " << u << ", " << v << " )" << std::endl;
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
                                                //std::cout << " T3 e: ( " << u << ", " << v << " )" << std::endl;
                                           count d_u = reduced_degree[u];
                                           count d_v = reduced_degree[v];
                                           total += sorted_intersection(u, edges2[u], d_u, v, edges2[v], d_v, m);
                                   });
                return total;
        }



                


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
