// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"
#include "math.h"

typedef double fType;

template <class vertex>
struct PR_F {
  fType* p_curr, *p_next;
  fType** p_next_dup;
  vertex* V;
  Bitmap& hubMap;
  uintE* idxMap;
  uintE* inv_idxMap;
  uintE hubCount;

  PR_F(fType* _p_curr, fType* _p_next, fType** _p_next_dup, vertex* _V, Bitmap& _hubMap, uintE* _idxMap, uintE* _inv_idxMap, uintE _hubCount ) : 
    p_curr(_p_curr), p_next(_p_next), p_next_dup(_p_next_dup), V(_V), hubMap(_hubMap), idxMap(_idxMap), inv_idxMap(_inv_idxMap), hubCount(_hubCount) {}
  inline bool update(uintE s, uintE d){ //update function applies PageRank equation
    p_next[d] += p_curr[s]/V[s].getOutDegree();
    return 0;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    if (hubMap.get_bit(d))
    {
        int tid = omp_get_thread_num();
        p_next_dup[tid][idxMap[d]] += (p_curr[s] / V[s].getOutDegree());
    }
    else
    {
        writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
    }
    return false;
  }

  /* merges thread-private duplicates. Also, resets thread-private 
  copies */
  inline void mergeDuplicates()
  {
      int numThreads = omp_get_max_threads();
      #pragma omp parallel for //schedule(static, 64/sizeof(fType))  //remove false sharing 
      for (intT d = 0; d < hubCount; ++d)
      {
          for (int tid = 0; tid < numThreads; ++tid) 
          {
              p_next[inv_idxMap[d]] += p_next_dup[tid][d];
              p_next_dup[tid][d] = 0.0;
          }
      }
      return;
  }

  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  fType damping;
  fType addedConstant;
  fType* p_curr;
  fType* p_next;
  PR_Vertex_F(fType* _p_curr, fType* _p_next, fType _damping, intE n) :
    p_curr(_p_curr), p_next(_p_next), 
    damping(_damping), addedConstant((1-_damping)*(1/(fType)n)){}
  inline bool operator () (uintE i) {
    p_next[i] = damping*p_next[i] + addedConstant;
    return 1;
  }
};

//resets p
struct PR_Vertex_Reset {
  fType* p_curr;
  PR_Vertex_Reset(fType* _p_curr) :
    p_curr(_p_curr) {}
  inline bool operator () (uintE i) {
    p_curr[i] = 0.0;
    return 1;
  }
};

void writeOutputToFile(fType* arr, long numElements, pvector<uintE> &new_ids) {
    std::ofstream fout;
    std::string prefixStr {"AppOutput.out"};
    fout.open(prefixStr);
    if (new_ids[0] == new_ids[1]) {
        // graph not preprocessed
        for (long v = 0; v < numElements; ++v) {
            fout << arr[v] << "\n";
        }
    }
    else {
        for (long v = 0; v < numElements; ++v) {
            fout << arr[new_ids[v]] << "\n";
        }
    }
    fout.close();
    return;
}


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, pvector<uintE> &new_ids) {
  Timer tm;
  tm.Start();
  long maxIters = P.getOptionLongValue("-maxiters",100);
  const intE n = GA.n;
  const fType damping = 0.85, epsilon = 0.0000001;
  //const fType damping = 0.85, epsilon = 0.0001;

  
  fType one_over_n = 1/(fType)n;
  #ifndef ALIGNED
  fType* p_curr = newA(fType,n);
  #else
  fType* p_curr {nullptr};
  posix_memalign((void**) &p_curr, 64, sizeof(fType) * n);
  assert(p_curr != nullptr && ((uintptr_t)p_curr % 64 == 0) && "App Malloc Failure\n");
  #endif
  {parallel_for(long i=0;i<n;i++) p_curr[i] = one_over_n;}
  fType* p_next = newA(fType,n);
  {parallel_for(long i=0;i<n;i++) p_next[i] = 0;} //0 if unchanged
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}

  vertexSubset Frontier(n,n,frontier);

  uintT avgDegree = GA.m / GA.n;
  //intT numHubs    = findNumHubs(GA.V, n, avgDegree); 
  intT llcCap   = 1 << 25; //32MB conservative so that all duplicates of hubs fit in our 35MB LLC
  size_t addCost = (sizeof(uintE) + sizeof(bool));  //bitmap cost should be 1/8B per hub (best-case)
                                                    //or 8B per hub (worst-case). Assuming GMean cost of 1B
  intT numHubs  = llcCap / ((omp_get_max_threads() * sizeof(fType)) + addCost);
  //create per-thread duplicates only for the hubvertices;
  int numThreads  = omp_get_max_threads();
  fType** p_next_dup = new fType* [numThreads];
  assert(p_next_dup != nullptr);
  for (int tid = 0; tid < numThreads; ++tid) 
  {
    p_next_dup[tid] = new fType [numHubs]();
    assert(p_next_dup[tid] != nullptr);
  }

  long iter = 0;
  fType L1_norm {0.0};
  
  //Populating data structures for efficient merge 
  tm.Stop();
  uintE* idxMap {nullptr};
  uintE* inv_idxMap {nullptr};
  Bitmap hubMap(GA.n);
  hubMap.reset();
  uintE hubCount {0};
  identifyHubs(GA, hubMap, idxMap, inv_idxMap, sizeof(fType), hubCount, true); 
  tm.Start();


  while(iter++ < maxIters) {
    bool* nextBitmap {nullptr};
    edgeMap(GA,Frontier,nextBitmap,PR_F<vertex>(p_curr,p_next,p_next_dup,GA.V,hubMap,idxMap,inv_idxMap,hubCount),-1, no_output | dense_forward);
    vertexMap(Frontier,PR_Vertex_F(p_curr,p_next,damping,n));
    //compute L1-norm between p_curr and p_next
    {parallel_for(long i=0;i<n;i++) {
      p_curr[i] = fabs(p_curr[i]-p_next[i]);
      }}
    L1_norm = sequence::plusReduce(p_curr,n);
    if(L1_norm < epsilon) break;
    //reset p_curr
    vertexMap(Frontier,PR_Vertex_Reset(p_curr));
    swap(p_curr,p_next);
  }
  //std::cout << "Num Iters until convergence = " << iter << std::endl;
  //writeOutputToFile(p_next, n, new_ids);
  Frontier.del(); free(p_curr); free(p_next); 
  for (int tid = 0; tid < numThreads; ++tid) 
  {
    delete[] p_next_dup[tid];
  }
  delete[] p_next_dup;
  tm.Stop();
  tm.PrintTime("Run Time(sec) ", tm.Seconds()); 
  std::cout << "[OUTPUT] Num Iters = " << iter << std::endl;
  std::cout << "[OUTPUT] L1_Norm   = " << L1_norm << std::endl;
  delete[] idxMap;
  delete[] inv_idxMap;
}
