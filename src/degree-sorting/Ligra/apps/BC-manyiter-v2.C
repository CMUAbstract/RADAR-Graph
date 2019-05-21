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
#include <vector>
#include <random>

typedef double fType;

const int NUM_ITERS {16};

struct BC_F {
  bool* nextBitmap;
  fType* NumPaths;
  bool* Visited;
  BC_F(bool* _nextBitmap, fType* _NumPaths, bool* _Visited) : 
    nextBitmap(_nextBitmap), NumPaths(_NumPaths), Visited(_Visited) {}
  inline bool update(uintE s, uintE d){ //Update function for forward phase
    fType oldV = NumPaths[d];
    NumPaths[d] += NumPaths[s];
    if (oldV == 0.0)
        nextBitmap[d] = 1;
    return false;
  }

  inline bool updateAtomic (uintE s, uintE d) { //atomic Update, basically an add
    volatile fType oldV, newV; 
    do { 
      oldV = NumPaths[d]; newV = oldV + NumPaths[s];
    } while(!CAS(&NumPaths[d],oldV,newV));
    if (oldV == 0.0)
        nextBitmap[d] = 1;
    return false;
  }
  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

struct BC_Back_F {
  fType* Dependencies;
  bool* Visited;
  BC_Back_F(fType* _Dependencies, bool* _Visited) : 
    Dependencies(_Dependencies), Visited(_Visited) {}
  inline bool update(uintE s, uintE d){ //Update function for backwards phase
    fType oldV = Dependencies[d];
    Dependencies[d] += Dependencies[s];
    return false;
  }

  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    volatile fType oldV, newV;
    do {
      oldV = Dependencies[d];
      newV = oldV + Dependencies[s];
    } while(!CAS(&Dependencies[d],oldV,newV));
    return false; 
  }

  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

//vertex map function to mark visited vertexSubset
struct BC_Vertex_F {
  bool* Visited;
  BC_Vertex_F(bool* _Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

//vertex map function (used on backwards phase) to mark visited vertexSubset
//and add to Dependencies score
struct BC_Back_Vertex_F {
  bool* Visited;
  fType* Dependencies, *inverseNumPaths;
  BC_Back_Vertex_F(bool* _Visited, fType* _Dependencies, fType* _inverseNumPaths) : 
    Visited(_Visited), Dependencies(_Dependencies), inverseNumPaths(_inverseNumPaths) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    Dependencies[i] += inverseNumPaths[i];
    return 1; }};

void writeOutputToFile(double* arr, long numElements, pvector<uintE> &new_ids) {
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
  /* App starts now */
  Timer tm;
  tm.Start();
  ulong sourceCtr {0};
  int outputIters {0};
  int maxIters = P.getOptionIntValue("-maxiters", NUM_ITERS);
  long n = GA.n;
  
  /* Allocating main data structures */
  #ifndef ALIGNED
  fType* NumPaths = newA(fType,n);
  fType* Dependencies = newA(fType,n);
  #else
  fType* NumPaths {nullptr};
  posix_memalign((void**) &NumPaths, 64, sizeof(fType) * n);
  assert(NumPaths != nullptr && ((uintptr_t)NumPaths % 64 == 0) && "App Malloc Failure\n");
  double* Dependencies {nullptr};
  posix_memalign((void**) &Dependencies, 64, sizeof(double) * n);
  assert(Dependencies != nullptr && ((uintptr_t)Dependencies % 64 == 0) && "App Malloc Failure\n");
  #endif
  {parallel_for(long i=0;i<n;i++) Dependencies[i] = 0.0;}
    
  bool* Visited = newA(bool,n);

  /* set up the random number generator */
  std::mt19937 rng(27491095);
  std::uniform_int_distribution<long> udist(0,GA.n-1);

  // Start the main computation - BC over a few samples
  for (int iteration = 0; iteration < maxIters; ++iteration) {
    bool preprocessed = (new_ids[0] != new_ids[1]); 
    long start = P.getOptionLongValue("-r",0);
    if (preprocessed) {
        do {
            //start = hashInt(sourceCtr++) % n; //random numbers are in original vertex space
            start = udist(rng);
        } while (GA.V[new_ids[start]].getOutDegree() == 0);
    }
    else {
        do {
            //start = hashInt(sourceCtr++) % n;
            start = udist(rng);
        } while (GA.V[start].getOutDegree() == 0);
    }

    #if 0
    tm.Stop();
    std::cout << "Source = " << start << std::endl;
    tm.Start();
    #endif
    
    if (preprocessed) 
      start = new_ids[start];

    
    {parallel_for(long i=0;i<n;i++) NumPaths[i] = 0.0;}
    NumPaths[start] = 1.0;
    {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
    Visited[start] = 1;
    vertexSubset Frontier(n,start);
 
    vector<vertexSubset> Levels;
    Levels.push_back(Frontier);

    long round = 0;
    while(!Frontier.isEmpty()){ //first phase
      round++;
      bool *nextBitmap = newA(bool, n);
      assert(nextBitmap != nullptr);
      vertexSubset output = edgeMap(GA, Frontier, nextBitmap, BC_F(nextBitmap,NumPaths,Visited), -1, dense_forward);
      vertexMap(output, BC_Vertex_F(Visited)); //mark visited
      Levels.push_back(output); //save frontier onto Levels
      Frontier = output;
    }
    outputIters += round;
    #if 0
    tm.Stop();
    std::cout << "[OUTPUT] First phase ends with " << round << " iterations" << std::endl;
    tm.Start();
    #endif

    //invert numpaths
    fType* inverseNumPaths = NumPaths;
    {parallel_for(long i=0;i<n;i++) inverseNumPaths[i] = 1/inverseNumPaths[i];}

    Levels[round].del();
    //reuse Visited
    {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
    Frontier = Levels[round-1];
    vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

    //tranpose graph
    GA.transpose();
    for(long r=round-2;r>=0;r--) { //backwards phase
      bool* nextBitmap {nullptr};
      edgeMap(GA, Frontier, nextBitmap, BC_Back_F(Dependencies,Visited), -1, no_output | dense_forward);
      Frontier.del();
      Frontier = Levels[r]; //gets frontier from Levels array
      //vertex map to mark visited and update Dependencies scores
      vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
    }
    
    Frontier.del();

    //Update dependencies scores
    parallel_for(long i=0;i<n;i++) {
      Dependencies[i]=(Dependencies[i]-inverseNumPaths[i])/inverseNumPaths[i];
    }
    //writeOutputToFile(Dependencies, n, new_ids);
  }
  free(NumPaths);
  free(Visited);
  free(Dependencies);
  tm.Stop();
  tm.PrintTime("Run Time(sec) ", tm.Seconds()); 
  std::cout << "[OUTPUT] Combined Iters (in front traversal) = " << outputIters << std::endl;
}
