/*
 * TriangleCounting.hpp
 *
 *  Created on: 29.04.21
 */

#ifndef NETWORKIT_GLOBAL_TRIANGLE_COUNTING_HPP_
#define NETWORKIT_GLOBAL_TRIANGLE_COUNTING_HPP_

#include <omp.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * @ingroup global
 */

class TriangleCounting : public Algorithm {
public:
        /**
         * Creates an instance of TriangleCounting for the given Graph @a graph.
         * @param graph
         */
        TriangleCounting(Graph &G) : G(&G) {
                // TODO: check implementation for directed graphs also.
                if (G.isDirected())
                        throw std::runtime_error("Error, Triangle Counting is not for directed graphs.");
                //triangles.assign(G->upperNodeIdBound(), none);
                total_count = 0;
        }
        ~TriangleCounting() override = default;

        /**
         * Computes the total number of triangles.
         */
        void run() override;
        void run_par();
        void run_seq();
        
        count getTriangleCount() {
                assureFinished();
                return total_count;
        }

        std::vector<count> getTriangleScores() {
                assureFinished();
                return triangles;
        }

        

private:
        // Graph& G;
        Graph *G; // TODO: const?
        std::vector<count> triangles;
        count total_count;


        void workPerThread(const count uLength, const count vLength, 
                           const int threadsPerIntersection, const int threadId,
                           int &outWorkPerThread, int & outDiagonalId);
        
        void initialize(const int diag_id, const count u_len, count v_len,
                        node &u_min, node &u_max,
                        node &v_min, node &v_max,
                        unsigned int &found );
        
        /* Requirement: u_len < v_len */
        count s_count_triangles(node u, const std::vector<node> & u_nodes, count u_len,
                                node v, const std::vector<node> & v_nodes, count v_len,
                                volatile node *  firstFound,
                                int tId
                                );
        
        count s_count_triangles_fg(node u, const std::vector<node> & u_nodes, count u_len,
                                   node v, const std::vector<node> & v_nodes, count v_len,
                                   volatile node *  firstFound,
                                   int tId
                                   );
        
        void s_intersect(const count uLength, const count vLength,
                         const std::vector<node> & uNodes, const std::vector<node> & vNodes,
                         node uCurr, node vCurr,
                         int  workIndex, int  workPerThread,
                         count &triangles, int found);
        

        void binary_search(unsigned int found, const int diagonalId,
                           const std::vector<node> & uNodes,   const std::vector<node> & vNodes,
                           count uLength, 
                           node outUMin, node outUMax,
                           node outVMin, node outVMax,    
                           node & outUCurr,
                           node & outVCurr);

        count intersections_seq (node u, node v, const std::vector<node> & uNodes, count uLength);
        
        
};

void TriangleCounting::run() {
        bool parallel = true;
        if (parallel) this->run_par();
        else run_seq();
}

count TriangleCounting::intersections_seq(node u, node v, const std::vector<node> & uNodes, count uLength) {

        count out = 0;
        G->forEdgesOf(v, [&](node, node w) {
                                 
                                 if (w != u) {
                                         int64_t first = 0;
                                         int64_t last = uLength-1;
                                         int64_t middle = (first + last)/2;
                                         
                                         while (first <= last) {
                                                 if (uNodes[middle] < w) {
                                                         first = middle + 1;
                                                 } else if (uNodes[middle] == w) {
                                                   out++;
                                                   break;
                                                 } else {
                                                         last = middle - 1;
                                                 }                 
                                                 middle = (first + last)/2;
                                         }
                                 }     
                         });
        return out;
}
        
void TriangleCounting::run_seq() {

        std::cout << " Running seq version ! " << std::endl;        
        std::fill(triangles.begin(), triangles.end(), 0);
        triangles.resize(G->upperNodeIdBound(), 0);
        
        std::vector<std::vector<node> > edges(G->upperNodeIdBound());
        /* sort edge lists */ 
        G->sortEdges();
        /* copy edge lists */ 
        G->parallelForNodes([&](node u) {
                                   edges[u].reserve(G->degree(u));
                                   G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                           edges[u].emplace_back(v);
                                                   });
                           });
        
        G->balancedParallelForNodes([&](node src) {                                            
                                            count srcLen = G->degree(src); 
                                            count tCount = 0;	    
                                            assert(!edges[src].empty());
                                            G->forEdgesOf(src, [&](node, node dest) {
                                                                       if ( dest != src ) {
                                                                               tCount += intersections_seq (src, dest,
                                                                                                            edges[src], srcLen);
                                                                       }
                                                               });
                                            triangles[src] = tCount;
                                    });
        total_count = G->parallelSumForNodes([&](node u) {
                                                     return triangles[u];
                                             });
        hasRun = true;
}


void TriangleCounting::run_par() {

        std::cout << " Running parallel version ! " << std::endl;

        std::fill(triangles.begin(), triangles.end(), 0);
        triangles.resize(G->upperNodeIdBound(), 0);
        
        std::vector<count> outPutTriangles(G->upperNodeIdBound(), 0);
        int tsp = 1; // Threads per intersection
        int sps = 128; // Block size
        // int nbl = sps/tsp; // Number of concurrent intersections in block
        // int blocks = 16000; // Number of blocks      
        const int threads_per_block = tsp;
        // const int number_blocks = nbl;
        const int shifter = 0; // left shift to multiply threads per intersection
        // const int thread_blocks = blocks;
        // const int blockdim = sps;
        // int tx = 0;
        // count this_mp_start, this_mp_stop;
        // this_mp_stop = n;
	// const int blockSize = blockDim.x;
	// workPerBlockIP(nv, &this_mp_start, &this_mp_stop, blockSize);

        std::vector<std::vector<node> > edges(G->upperNodeIdBound());
        G->sortEdges();
        G->parallelForNodes([&](node u) {
                                   edges[u].reserve(G->degree(u));
                                   G->forEdgesOf(u, [&](node, node v, edgeid) {
                                                           edges[u].emplace_back(v);
                                                   });
                           });

        //std::vector<std::vector<count> > incidentTriangleCount(omp_get_max_threads(), std::vector<count>(G->upperNodeIdBound(), none));
	count s_triangles[1024];
        node firstFound[1024];
        
        const count maxThreads = static_cast<count>(omp_get_max_threads());

        int adj_offset;         
        index tid;
        long unsigned int * firstFoundPos;

#pragma omp parallel private(tid, adj_offset, firstFoundPos)
        {
                // auto tid = omp_get_thread_num();
                tid = omp_get_thread_num();
                adj_offset = tid >> shifter;
                firstFoundPos = firstFound + (adj_offset << shifter);
                usleep(5000 * omp_get_thread_num()); // do this to avoid race condition while printing
                std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;
                // each thread can also get its own number
                std::cout << "Current thread number: " << omp_get_thread_num() << " adj_offset = " << adj_offset << std::endl;
                std::cout << "Hello, World!" << std::endl;
        }
        // for all nodes
        G->balancedParallelForNodes([&](node src) {
                                            count srcLen = G->degree(src); 
                                            count tCount = 0;	    
                                            //G->forEdgesOf(src, [&](node, node dest) {
                                            assert(!edges[src].empty());
                                            //std::cout << "source : " << src << " d : " << srcLen << "\n edges: " << std::endl;
                                            //for (index i = 0; i < edges[src].size(); ++i) std::cout << " " << edges[src][i];
                                            //std::cout << std::endl;                            
                                            for (index i = 0; i < edges[src].size(); ++i) {                          
                                                    node dest = edges[src][i];
                                                    count destLen = G->degree(dest);
                                                    bool avoidCalc = (src == dest) || (destLen < 2) || (srcLen < 2);
                                                    if(avoidCalc)
                                                            continue;
                                    
                                                    bool sourceSmaller = (srcLen<destLen);
                                                    node small = sourceSmaller? src : dest;
                                                    node large = sourceSmaller? dest : src;
                                                    count small_len = sourceSmaller? srcLen : destLen;
                                                    count large_len = sourceSmaller? destLen : srcLen;
                                    
                                                    const std::vector<node> small_ptr = edges[small];
                                                    const std::vector<node> large_ptr = edges[large];
                                                    assert(!large_ptr.empty() && !small_ptr.empty() );
                                                    
                                                    tCount += s_count_triangles(small, small_ptr, small_len,
                                                                                large, large_ptr, large_len,
                                                                                //threads_per_block,
                                                                                firstFoundPos, tid%threads_per_block);
                                            }
                                    //});
                            // adjust for more than one threads
                            s_triangles[tid] = tCount;
                            //blockReduceIP(&outPutTriangles[src],s_triangles,blockSize);
                    });


        count t_t = 0;        
#pragma omp parallel for reduction(+ : t_t)
        for (omp_index v = 0; v < static_cast<omp_index>(maxThreads); ++v)
                if (s_triangles[v])
                        t_t += s_triangles[v];        
#pragma omp parallel for
        for (index i = 0; i < maxThreads; ++i)
                std::cout << " s_triangles of thread # " << tid << " = " << s_triangles[i] << std::endl;

        total_count = t_t;

        // for (node u = 0; u < triangles.size(); u++)
        //         total_count += triangles[u]; 

        // or:
        // double cc = G.parallelSumForNodes([&](node u){
        //                                           return triangles[u];
        //                                   });
        hasRun = true;

}
        

// void TriangleCounting::run() {
        
//         std::fill(triangles.begin(), values.end(), 0);
//         triangles.resize(G.upperNodeIdBound(), 0);
// //        cuStinger* custing;
// //        void callDeviceAllTriangles(cuStinger& custing,
// //                                    triangle_t * const __restrict__ outPutTriangles, const int threads_per_block,
// //                                    const int number_blocks, const int shifter, const int thread_blocks, const int blockdim){
// //        devicecuStingerAllTriangles<<<thread_blocks, blockdim>>>(custing.devicePtr(), outPutTriangles,
// //                                                                 threads_per_block,number_blocks,shifter); }

        
//         std::vector<count> outPutTriangles(G.upperNodeIdBound(), 0);
//         int tsp = 1; // Threads per intersection
//         int sps = 128; // Block size
//         int nbl = sps/tsp; // Number of concurrent intersections in block
//         int blocks = 16000; // Number of blocks
        
//         const int threads_per_block = tsp;
//         const int number_blocks = nbl;
//         const int shifter = 0; // left shift to multiply threads per intersection
//         const int thread_blocks = blocks;
//         const int blockdim = sps;

        
// 	count n = G.upperNodeIdBound();
// 	// Partitioning the work to the multiple threads.
//         // The threads should get a near equal number of the elements to intersect - this number will be off by no more than one.
// 	// int tx = threadIdx.x;
//         int tx = 0;
//         count this_mp_start, this_mp_stop;
//         this_mp_stop = n;

// 	// const int blockSize = blockDim.x;
// 	// workPerBlockIP(nv, &this_mp_start, &this_mp_stop, blockSize);
        
// 	count s_triangles[1024];
//         node firstFound[1024];

// 	int adj_offset=tx>>shifter;
//         int * firstFoundPos=firstFound + (adj_offset<<shifter);
//         // for all nodes
//         for (node src = this_mp_start; src < this_mp_stop; src++){
// 		// int srcLen=d_off[src+1]-d_off[src];
// 		count srcLen= G.degree(src); 
//                 count tCount = 0;	    
// 		// for(int iter=d_off[src]+adj_offset; iter<d_off[src+1]; iter+=number_blocks){
// 		for(int k=adj_offset; k<srcLen; k+=number_blocks){
// 			// // int dest = d_ind[k];
// 			// vertexId_t dest = custing->dVD->getAdj()[src]->dst[k];
//                         node dest = G.getTarget(src);
// 			// // int destLen = d_off[dest+1]-d_off[dest];
// 			// int destLen=custing->dVD->getUsed()[dest];
//                         count destLen = G.degree(dest);
                        
// 			// if (dest<src) 
// 			// 	continue;

// 			bool avoidCalc = (src == dest) || (destLen < 2) || (srcLen < 2);
// 			if(avoidCalc)
// 				continue;

//                         bool sourceSmaller = (srcLen<destLen);
//                         node small = sourceSmaller? src : dest;
//                         node large = sourceSmaller? dest : src;
//                         count small_len = sourceSmaller? srcLen : destLen;
//                         count large_len = sourceSmaller? destLen : srcLen;


// 	        // // int const * const small_ptr = d_ind + d_off[small];
// 	        // // int const * const large_ptr = d_ind + d_off[large];
// 	        // const vertexId_t* small_ptr = custing->dVD->getAdj()[small]->dst;
// 	        // const vertexId_t* large_ptr = custing->dVD->getAdj()[large]->dst;
                        
//                         const std::vector<node> small_ptr = G.getAdj(small);
//                         const std::vector<node> large_ptr = G.getAdj(large);
//                         assert(!large_ptr.empty() && !small_ptr.empty() );
                        
//                         tCount += s_count_triangles(small, small_ptr, small_len,
//                                                     large, large_ptr, large_len,
//                                                     //threads_per_block,
//                                                     firstFoundPos, tx%threads_per_block);
// 		}
// 		s_triangles[tx] = tCount;
// 		//blockReduceIP(&outPutTriangles[src],s_triangles,blockSize);
// 	}
//         for (node u = 0; u < triangles.size(); u++){
//                 total_count += triangles[u]; 
//         }
//         // or:
//         // double cc = G.parallelSumForNodes([&](node u){
//         //                                           return triangles[u];
//         //                                   });

//         hasRun = true;

// }
        


        /* count triangles */
        count TriangleCounting::s_count_triangles_fg(node u, const std::vector<node> & u_nodes, count u_len,
                                                     node v, const std::vector<node> & v_nodes, count v_len,
                                                     // int threads_per_block,
                                                     // TODO : std::vector<node>
                                                     volatile node *  firstFound,
                                                     int tId
                                                     ) {
                // Partitioning the work to the multiple thread of a single GPU processor.
                // The threads should get a near equal number of the elements (off by 1).  
                int work_per_thread, diag_id;
                const int threads_per_block = 1; 
                workPerThread(u_len, v_len, threads_per_block, tId, work_per_thread, diag_id);
                if (!diag_id) std::cout << "Diag id is not zero! Ok " << std::endl;
                if (work_per_thread != (u_len + v_len)) std::cout << " Work not the entire lenght! Weird " << std::endl;
                //assert( work_per_thread != (u_len + v_len) );
                
                count triangles = 0;
                int work_index = 0;
                unsigned int found = 0;
                node u_min, u_max, v_min, v_max, u_curr, v_curr;

                // TODO: what is firstFound ? Why is it volatile in the cuda code?
                firstFound[tId] = 0;
                std::cout << "pair (" << u << ", " << v << ")  -> work_per_thread = " << work_per_thread << std::endl;
                
                if(work_per_thread > 0){
                        // For the binary search, find the initial poT of search.
                        initialize(diag_id, u_len, v_len, u_min, u_max, v_min, v_max, found);
                        u_curr = 0; v_curr = 0;                        
                        binary_search(found, diag_id, u_nodes, v_nodes, u_len, u_min, u_max, v_min,
                                v_max, u_curr, v_curr);
                        
                        // fix starting point
                        node uBigger = (u_curr > 0) && (v_curr < v_len) && (u_nodes[u_curr-1] == v_nodes[v_curr]);
                        node vBigger = (v_curr > 0) && (u_curr < u_len) && (v_nodes[v_curr-1] == u_nodes[u_curr]);
                        u_curr += vBigger;
                        v_curr += uBigger;
                        int sum = uBigger + vBigger;
                        // end fix starting
                        work_index += sum;
                        
                        if(tId > 0)
                                firstFound[tId-1] = sum;
                        triangles += sum;
                      
                        s_intersect(u_len, v_len, u_nodes, v_nodes, u_curr, v_curr,
                                    work_index, work_per_thread, triangles, firstFound[tId]);
                }
                return triangles;
        }


        
       
        /* count triangles */
        count TriangleCounting::s_count_triangles(node u, const std::vector<node> & u_nodes, count u_len,
                                                  node v, const std::vector<node> & v_nodes, count v_len,
                                                  // int threads_per_block,
                                                  // TODO : std::vector<node>
                                                  volatile node *  firstFound,
                                                  int tId
                                                  ) {
                int work_per_thread = u_len + v_len;
                int diag_id = 0;
                // END TODO
                
                count triangles = 0;
                int work_index = 0;
                unsigned int found = 0;
                node u_min, u_max, v_min, v_max, u_curr, v_curr;

                // TODO: what is firstFound ? Why is it volatile in the cuda code?
                firstFound[tId] = 0;
                std::cout << "pair (" << u << ", " << v << ")  -> work_per_thread = " << work_per_thread << std::endl;
                
                if(work_per_thread > 0){
                        // 	// For the binary search, we are figuring out the initial poT of search.
                        // 	initialize(diag_id, u_len, v_len,&u_min, &u_max,&v_min, &v_max,&found);
                        // Currently for one thread -> TODO: change
                        u_min = u_max = v_min= v_max=0;
                        found = 1;
                        // END TODO
                        u_curr = 0; v_curr = 0;
                        
                        binary_search(found, diag_id, u_nodes, v_nodes, u_len, u_min, u_max, v_min,
                                v_max, u_curr, v_curr);
                        
                        // fix starting point
                        node uBigger = (u_curr > 0) && (v_curr < v_len) && (u_nodes[u_curr-1] == v_nodes[v_curr]);
                        node vBigger = (v_curr > 0) && (u_curr < u_len) && (v_nodes[v_curr-1] == u_nodes[u_curr]);
                        u_curr += vBigger;
                        v_curr += uBigger;
                        int sum = uBigger + vBigger;
                        // end fix starting
                        work_index += sum;
                        
                        if(tId > 0)
                                firstFound[tId-1] = sum;
                        triangles += sum;
                        
                        s_intersect(u_len, v_len, u_nodes, v_nodes, u_curr, v_curr,
                                    work_index, work_per_thread, triangles, firstFound[tId]);
                }
                return triangles;
        }
        
                 
        /* intersect adjacency lists and 'return' triangles */
        /* const std::vector<node> & uNodes, & vNodes better than just const std::vector<node> uNodes, vNodes ?  */
        void TriangleCounting::s_intersect(const count uLength, const count vLength,
                                           const std::vector<node> & uNodes, const std::vector<node> & vNodes,
                                           node uCurr, node vCurr,
                                           int  workIndex, int  workPerThread,
                                           count &triangles, int found) {
                
                assert(!uNodes.empty() && !vNodes.empty());
                assert( (uNodes.size() == uLength) && (vNodes.size() == vLength) );
                //std::cout << "cur (" << uCurr << ", " << vCurr << ")  -> len ( " << uLength << ", " << vLength << ") " << std::endl;
                //std::cout << " workIndex =  " << workIndex << ", workPerThread = " << workPerThread << std::endl;
                if( (uCurr < uLength) && (vCurr < vLength) ) {
                        int comp;                
                        while( workIndex < workPerThread ) {
                                comp = uNodes[uCurr] - vNodes[vCurr];                        
                                triangles  += (comp == 0);
                                
                                uCurr += (comp <= 0);
                                vCurr += (comp >= 0);
                                workIndex += (comp == 0) + 1;
                                
                                if((vCurr == vLength) || (uCurr == uLength)){
                                        break;
                                }
                        }
                        
                        // std::cout << " c = " <<  (comp == 0) << std::endl; 
                        // std::cout << " w = " <<  (workIndex > workPerThread) << std::endl; 
                        // std::cout << " f = " <<  found << std::endl;
                        // std::cout << " t = " <<  ((comp == 0) && (workIndex > workPerThread) && (found)) << std::endl;
                        triangles -= ((comp == 0) && (workIndex > workPerThread) && (found));
                }
        }

        /* perform binary search */
        void TriangleCounting::binary_search(unsigned int found, const int diagonalId,
                                             const std::vector<node> & uNodes, const std::vector<node> & vNodes,
                                             count uLength, 
                                             node outUMin, node outUMax,
                                             node outVMin, node outVMax,    
                                             node &outUCurr,
                                             node &outVCurr) {
                
                count length;
                
                while(!found) {
                        outUCurr = (outUMin + outUMax)>>1;
                        outVCurr = diagonalId - outUCurr;
                        if(outVCurr >= outVMax){
                                length = outUMax - outUMin;
                                if(length == 1){
                                        found = 1;
                                        continue;
                                }
                        }
                     
                        unsigned int comp1 = uNodes[outUCurr] > vNodes[outVCurr-1];
                        unsigned int comp2 = uNodes[outUCurr-1] > vNodes[outVCurr];
                        if(comp1 && !comp2){
                                found = 1;
                        }
                        else if(comp1){
                                outVMin = outVCurr;
                                outUMax = outUCurr;
                        }
                        else{
                                outVMax = outVCurr;
                                outUMin = outUCurr;
                        }
                }
                
                if((outVCurr >= outVMax) && (length == 1) && (outVCurr > 0) &&
                   (outUCurr > 0) && (outUCurr < (uLength - 1))){
                        unsigned int comp1 = uNodes[outUCurr] > vNodes[outVCurr - 1];
                        unsigned int comp2 = uNodes[outUCurr - 1] > vNodes[outVCurr];
                        if(!comp1 && !comp2){ outUCurr++; outVCurr--;}
                }
        }



void TriangleCounting::workPerThread(const count uLength, const count vLength, 
                   const int threadsPerIntersection, const int threadId,
                   int &outWorkPerThread, int & outDiagonalId){
        int totalWork = uLength + vLength;
        int remainderWork = totalWork%threadsPerIntersection; // 0 for threadsPerIntersection = 1
        int workPerThread = totalWork/threadsPerIntersection; // totalWork  for tpi = 1
        
        int longDiagonals  = (threadId > remainderWork) ? remainderWork:threadId; // 0 for tpi = 1
        int shortDiagonals = (threadId > remainderWork) ? (threadId - remainderWork):0; // equal to threadId for tpi = 1
        
        outDiagonalId = ((workPerThread+1)*longDiagonals) + (workPerThread*shortDiagonals); // equal to workPerThread*threadId for tpi = 1
        outWorkPerThread = workPerThread + (threadId < remainderWork); // workPerThred for tpi = 1
}

void TriangleCounting::initialize(const int diag_id, const count u_len, count v_len,
                  node &u_min, node &u_max,
                  node &v_min, node &v_max,
                  unsigned int &found )
{
	if (diag_id == 0){
		u_min = u_max = v_min = v_max = 0;
		found=1;
	}
	else if (diag_id < u_len){
		u_min = 0; u_max = diag_id;
		v_max = diag_id; v_min = 0;
	}
	else if (diag_id < v_len){
		u_min=0; u_max = u_len;
		v_max=diag_id; v_min = diag_id-u_len;
	}
	else{
		u_min = diag_id-v_len; u_max = u_len;
		v_min = diag_id-u_len; v_max = v_len;
	}
}
        

} /* namespace NetworKit */


#endif // NETWORKIT_GLOBAL_TRIANGLE_COUNTING_HPP_
