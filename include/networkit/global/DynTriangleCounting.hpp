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
        
        count intersectionCount( node u, std::vector<node> & u_adj, bool u_red,
                                 node v, std::vector<node> & v_adj, bool v_red, double m);
        
        // consider call by ref for adj_inG and adj_inB
        count do_merge_insert(const std::vector<node> adj_inG, const count deg_inG,
                              const std::vector<node> adj_inB, const count deg_inB,
                              std::vector<node> & adj_inGI);

        count do_merge_delete(std::vector<node> & adj_inG, const count deg_inG,
                              const std::vector<node> & adj_inB, const count deg_inB);

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
                G->balancedParallelForNodes([&](node u) {
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
                // storing previous total count (due to static run)
                o_t = t_t;
                // check batch input, create batch update graph.
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
                // sort update graph.
                ugraph.sortEdges();

                // create local adjacency list copies for G
                // adjacency lists are used for merging graphs directly and counting triangles
                std::vector<std::vector<node> > edges(G->upperNodeIdBound());
                G->balancedParallelForNodes([&](node u) {
                                                    edges[u].reserve(G->degree(u));
                                                    G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                                             edges[u].emplace_back(v);
                                                                     });
                                            });
                
                // create local adjacency list copies for update graph
                std::vector<std::vector<node> > edges2(ugraph.upperNodeIdBound());
                // TODO: consider going through edges instead of nodes since m << n in update graph
                // TODO: if so, consider costs by locking simultaneous access to nodes 
                ugraph.balancedParallelForNodes([&](node u) {
                                                        if (!ugraph.degree(u)) return;
                                                        edges2[u].reserve(ugraph.degree(u));
                                                        ugraph.forEdgesOf(u, [&](node, node v, edgeid) {
                                                                                     edges2[u].emplace_back(v);
                                                                             });
                                                });

                

                Aux::Timer timer;
                if (graphAlreadyUpdated) INFO("** ***** GRAPH ALREADY UPDATED *****");
                else {
                        // code to merge two sorted graphs, orignal graph G and update graph 
                        if(!checkSorted(NULL)) {
                                WARN("** ***** GRAPH WAS NOT SORTED. SORTING NOW! *****");
                                G->sortEdges();
                        }
                        std::vector<edgeweight> empty;
                        
                        timer.start();
                        
                        GraphBuilder result(G->upperNodeIdBound(), false, false);
                        G->balancedParallelForNodes([&](node v) {
                                                            count current_degree = G->degree(v);
                                                            count update_degree = ugraph.degree(v);
                                                            count max_new_degree = current_degree + update_degree;
                                                            // isolated nodes in both graphs
                                                            if (!current_degree && !update_degree)
                                                                    return;                                                     
                                                            if (insertion) {
                                                                    // inserting edges
                                                                    if (update_degree && current_degree) {
                                                                            std::vector<node> new_edges(max_new_degree); // reserve ?
                                                                            count new_degree = do_merge_insert(edges[v], current_degree, edges2[v], update_degree, new_edges);
                                                                            assert(new_edges.size() == new_degree);
                                                                            // copy merged adjacency list
                                                                            edges[v] = new_edges;
                                                                            // TODO: current implementation for unweighted graphs
                                                                            result.assignNeighborhood(v, new_edges, empty, false);
                                                                    }
                                                                    else if (!update_degree) {
                                                                            // no changes to adjacency list
                                                                            result.assignNeighborhood(v, edges[v], empty, false);
                                                                    }
                                                                    else {
                                                                            // isolated node in G
                                                                            result.assignNeighborhood(v, edges2[v], empty, false);
                                                                            edges[v] = edges2[v];
                                                                            
                                                                    }
                                                                    
                                                            }
                                                            else {
                                                                    if (update_degree) {
                                                                            // deleting edges
                                                                            // do_merge_delete changes adjacency lists in place
                                                                            count new_degree = do_merge_delete(edges[v], current_degree, edges2[v], update_degree);
                                                                            assert(new_degree == edges[v].size());
                                                                    }
                                                                    result.assignNeighborhood(v, edges[v], empty, false);
                                                            }
                                                    });
                        
                        Graph GU = result.toGraph(false,true);
                        U = std::move(GU);
                        G = &U;

                        timer.stop();
                        INFO(" ** ***** Time to insert batch =  ", timer.elapsedMicroseconds() / 1e6, " secs.");  
                }
                
                // setting parameters for computing triangles via intersection
                double mtype1, mtype2, mtype3;
                if (insertion) {
                        // for inserting edges
                        mtype1 = 1.0;
                        mtype2 = -1.0;
                        mtype3 = 1.0;
                }
                else {
                        // for deleting edges
                        mtype1 = mtype2 = mtype3 = -1.0;
                        if (!t_t) return;
                }

                ugraph.parallelForEdges([&](node u, node v) {

                                                // work on updated graph only
                                                intersectionCount(u, edges[u], false, v, edges[v], false, mtype1); 
                                                // work on updated (reduced adjacency lists) and update graph (full)
                                                intersectionCount(u, edges[u], true, v, edges2[v], false, mtype2);
                                                intersectionCount(v, edges[v], true, u, edges2[u], false, mtype2);
                                                // // work on update graph only (reduced adjacency lists)
                                                intersectionCount(u, edges2[u], true, v, edges2[v], true, mtype3); 
                                        });
                
                computeTriangleCount();
        }

       

        count DynTriangleCounting::intersectionCount( node u, std::vector<node> & u_adj, bool u_red,
                                 node v, std::vector<node> & v_adj, bool v_red, double m) {

                if ( u_adj.empty() || v_adj.empty() ) return 0;
                count triangles = 0;
                node i = 0, j = 0;
                count u_deg, v_deg;
                node * u_ptr;
                node * v_ptr;
                (u_red) ? u_deg = u : u_deg = u_adj.size();
                (v_red) ? v_deg = v : v_deg = v_adj.size();
                (u_red) ? u_ptr = &u_adj[i] : u_ptr = &i;
                (v_red) ? v_ptr = &v_adj[j] : v_ptr = &j; 


                // std::cout << " e: ( " << u << ", " << v << " ) --> list ( [";
                // for (int k = 0; k < u_deg; k++)
                //         std::cout << u_adj[k] << " ";
                // std::cout << "] , [";
                // for (int k = 0; k < v_deg; k++)
                //         std::cout << v_adj[k] << " ";
                // std::cout << "]) "<<  std::endl;
                
                while( ( (*u_ptr < u_deg) && (i < u_adj.size()) ) && ( (*v_ptr < v_deg) &&  (j < v_adj.size()) ) ) {
                        int comp;
                        comp = u_adj[i] - v_adj[j];
                        triangles += (comp == 0);
                        if (comp == 0) {
#pragma omp atomic
                                TrianglesPerNode[u_adj[i]] += m;
                        }
                        i += (comp <= 0);
                        j += (comp >= 0);
                        if(u_red) u_ptr = &u_adj[i];
                        if(v_red) v_ptr = &v_adj[j]; 
                        // std::cout << " u_ptr = " << *u_ptr << "\t i = " << i << std::endl;
                        // make sure i,j dont surpass limits of u_adj and v_adj -> seg fault 
                        // if( (i >= u_adj.size() && u_red) || (j >= v_adj.size() && v_red)) {
                        //         break;
                        // }
                }
#pragma omp atomic
                TrianglesPerNode[u] += m * triangles;
#pragma omp atomic
                TrianglesPerNode[v] += m * triangles;
                return triangles;
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
                        // TODO : break
                        G->balancedParallelForNodes([&] (node u) {
                                                            if(!std::is_sorted(G->neighborRange(u).begin(), G->neighborRange(u).end()))
                                                                    res = false;
                                                    });
                }
                else {
                        G1->balancedParallelForNodes([&] (node u) {
                                                             if(!std::is_sorted(G1->neighborRange(u).begin(), G1->neighborRange(u).end()))
                                                                     res = false;
                                                     });
                }
                return res;
        }
                                   
                        



        
} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_DYN_TRIANGLE_COUNTING_HPP_
