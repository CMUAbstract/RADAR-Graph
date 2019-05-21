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
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <parallel/algorithm>
#include <omp.h>
#include <cassert>

#include "parallel.h"
#include "blockRadixSort.h"
#include "quickSort.h"
#include "utils.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"
#include "sliding_queue.h"
#include "bitmap.h"

using namespace std;

typedef pair<uintE,uintE> intPair;
typedef pair<uintE, pair<uintE,intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator() (pair<uintE,E> a, pair<uintE,E> b) {
    return a.first < b.first; }
};

template <class E>
struct getFirst {uintE operator() (pair<uintE,E> a) {return a.first;} };

template <class IntType>
struct pairBothCmp {
  bool operator() (pair<uintE,IntType> a, pair<uintE,IntType> b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n; // total number of characters
  char* Chars;  // array storing all strings
  long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long nn, char** S, long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};

inline bool isSpace(char c) {
  switch (c)  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

_seq<char> mmapStringFromFile(const char *filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG (sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char *p = static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
//  char *bytes = newA(char, n);
//  parallel_for(size_t i=0; i<n; i++) {
//    bytes[i] = p[i];
//  }
//  if (munmap(p, sb.st_size) == -1) {
//    perror("munmap");
//    exit(-1);
//  }
//  cout << "mmapped" << endl;
//  free(bytes);
//  exit(0);
  return _seq<char>(p, n);
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file (fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg (0, ios::beg);
  long n = end - file.tellg();
  char* bytes = newA(char,n+1);
  assert(bytes != NULL && "Malloc failure\n");
  file.read (bytes,n);
  file.close();
  return _seq<char>(bytes,n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long n) {
  {parallel_for (long i=0; i < n; i++)
      if (isSpace(Str[i])) Str[i] = 0; }

  // mark start of words
  bool *FL = newA(bool,n);
  assert(FL != NULL && "Malloc failure\n");
  FL[0] = Str[0];
  {parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];}

  // offset for each start of word
  _seq<long> Off = sequence::packIndex<long>(FL, n);
  free(FL);
  long m = Off.n;
  long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char*, m);
  assert(SA != NULL && "Malloc failure\n");
  {parallel_for (long j=0; j < m; j++) SA[j] = Str+offsets[j];}

  free(offsets); 
  return words(Str,n,SA,m);
}

/* prefix sum used by the preprocess function defined below */
static pvector<uintT> ParallelPrefixSum (const pvector<uintT> &degrees) {
  const size_t block_size = 1<<20;
  const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
  pvector<uintT> local_sums(num_blocks);
  #pragma omp parallel for
  for (size_t block=0; block < num_blocks; block++) {
    uintT lsum = 0;
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i=block * block_size; i < block_end; i++)
      lsum += degrees[i];
    local_sums[block] = lsum;
  }
  pvector<uintT> bulk_prefix(num_blocks+1);
  uintT total = 0;
  for (size_t block=0; block < num_blocks; block++) {
    bulk_prefix[block] = total;
    total += local_sums[block];
  }
  bulk_prefix[num_blocks] = total;
  pvector<uintT> prefix(degrees.size() + 1);
  #pragma omp parallel for
  for (size_t block=0; block < num_blocks; block++) {
    uintT local_total = bulk_prefix[block];
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i=block * block_size; i < block_end; i++) {
      prefix[i] = local_total;
      local_total += degrees[i];
    }
  }
  prefix[degrees.size()] = bulk_prefix[num_blocks];
  return prefix;
}

template <class vertex>
graph<vertex> readGraphFromFile(char* fname, bool isSymmetric, bool mmap) {
  Timer t;
  t.Start();
  words W;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    assert(bytes != NULL && "Malloc failure\n");
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph") {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2) {
#else
  if (len != n + 2*m + 2) {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  uintT* offsets = newA(uintT,n);
  assert(offsets != NULL && "Malloc failure\n");
#ifndef WEIGHTED
  uintE* edges = newA(uintE,m);
#else
  intE* edges = newA(intE,2*m);
#endif
  assert(edges != NULL && "Malloc failure\n");

  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
#ifndef WEIGHTED
      edges[i] = atol(W.Strings[i+n+3]);
#else
      edges[2*i] = atol(W.Strings[i+n+3]);
      edges[2*i+1] = atol(W.Strings[i+n+m+3]);
#endif
    }}
    
  //W.del(); // to deal with performance bug in malloc
  W.del(); //The original code ^ commented this out

  vertex* v = newA(vertex,n);
  assert(v != NULL && "Malloc failure\n");

  {parallel_for (uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
#ifndef WEIGHTED
    v[i].setOutNeighbors(edges+o);
#else
    v[i].setOutNeighbors(edges+2*o);
#endif
  }}

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    assert(tOffsets != NULL && "Malloc failure\n");
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m);
#else
    intTriple* temp = newA(intTriple,m);
#endif
    assert(temp != NULL && "Malloc failure\n");
    {parallel_for(long i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    assert(inEdges != NULL && "Malloc failure\n");
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);

    {parallel_for(long i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors(inEdges+o);
#else
      v[i].setInNeighbors(inEdges+2*o);
#endif
      }}

    free(tOffsets);
    
    #ifndef WEIGHTED
    /* Removing redundant edges and self-loops*/ 
    
    pvector<uintT> new_out_degrees(n);
    pvector<uintT> new_in_degrees(n);
    uintE *n_out_start, *n_out_end; 
    uintE *n_in_start, *n_in_end; 
    #pragma omp parallel for private(n_out_start, n_out_end, n_in_start, n_in_end)
    for (long i = 0; i < n; ++i)
    {
        n_out_start = v[i].getOutNeighbors();
        n_out_end   = v[i+1].getOutNeighbors();
        n_in_start  = v[i].getInNeighbors();
        n_in_end    = v[i+1].getInNeighbors();
        if (i == n-1)
        {
            n_out_end = v[0].getOutNeighbors() + m;
            n_in_end  = v[0].getInNeighbors() + m;
        }

        std::sort(n_out_start, n_out_end);
        std::sort(n_in_start, n_in_end);

        uintE* new_end = std::unique(n_out_start, n_out_end);
        new_end = std::remove(n_out_start, new_end, i);
        new_out_degrees[i] = new_end - n_out_start;

        new_end = std::unique(n_in_start, n_in_end);
        new_end = std::remove(n_in_start, new_end, i);
        new_in_degrees[i] = new_end - n_in_start;
    }

    pvector<uintT> offsets     = ParallelPrefixSum(new_out_degrees);
    pvector<uintT> inv_offsets = ParallelPrefixSum(new_in_degrees);
    
    long new_numEdges {0};
    #pragma omp parallel for reduction(+ : new_numEdges)
    for (long i = 0; i < n; ++i)
        new_numEdges += new_out_degrees[i];
    
    vertex* newV        = newA(vertex, n);
    uintE* new_outEdges = newA(uintE, new_numEdges);
    uintE* new_inEdges  = newA(uintE, new_numEdges);

    #pragma omp parallel for schedule(dynamic, 64)
    for (long i = 0; i < n; ++i)
    {
        newV[i].setOutDegree(new_out_degrees[i]);
        newV[i].setOutNeighbors(new_outEdges + offsets[i]);
        for (uintE ngh = 0; ngh < new_out_degrees[i]; ++ngh)
            newV[i].setOutNeighbor(ngh, v[i].getOutNeighbor(ngh));
        
        newV[i].setInDegree(new_in_degrees[i]);
        newV[i].setInNeighbors(new_inEdges + inv_offsets[i]);
        for (uintE ngh = 0; ngh < new_in_degrees[i]; ++ngh)
            newV[i].setInNeighbor(ngh, v[i].getInNeighbor(ngh));
    }
    #endif
         
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,n,new_numEdges,new_outEdges,new_inEdges);
    t.Stop();
    t.PrintTime("Graph reading time(s)", t.Seconds());
    std::cout << "Read directed graph. Num Nodes = " << n << " and Num Edges = " << new_numEdges << "\n";
    return graph<vertex>(newV,n,new_numEdges,mem);
  }
  else {
    free(offsets);
    
    #ifndef WEIGHTED
    /* Removing redundant edges - out edges first*/ 
    
    pvector<uintT> new_out_degrees(n);
    uintE *n_out_start, *n_out_end; 
    #pragma omp parallel for private(n_out_start, n_out_end)
    for (long i = 0; i < n; ++i)
    {
        n_out_start = v[i].getOutNeighbors();
        n_out_end   = v[i+1].getOutNeighbors();
        if (i == n-1) n_out_end = v[0].getOutNeighbors() + m;
        std::sort(n_out_start, n_out_end);
        uintE* new_end = std::unique(n_out_start, n_out_end);
        new_end = std::remove(n_out_start, new_end, i);
        new_out_degrees[i] = new_end - n_out_start;
    }

    pvector<uintT> offsets     = ParallelPrefixSum(new_out_degrees);
    
    long new_numEdges {0};
    #pragma omp parallel for reduction(+ : new_numEdges)
    for (long i = 0; i < n; ++i)
        new_numEdges += new_out_degrees[i];
    
    vertex* newV        = newA(vertex, n);
    uintE* new_outEdges = newA(uintE, new_numEdges);

    #pragma omp parallel for schedule(dynamic, 64)
    for (long i = 0; i < n; ++i)
    {
        newV[i].setOutDegree(new_out_degrees[i]);
        newV[i].setOutNeighbors(new_outEdges + offsets[i]);
        for (uintE ngh = 0; ngh < new_out_degrees[i]; ++ngh)
            newV[i].setOutNeighbor(ngh, v[i].getOutNeighbor(ngh));
    }
    #endif

    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,n,new_numEdges,new_outEdges);
    t.Stop();
    t.PrintTime("Graph reading time(s)", t.Seconds());
    std::cout << "Read undirected graph. Num Nodes = " << n << " and Num Edges = " << new_numEdges << "\n";
    return graph<vertex>(newV,n,new_numEdges,mem);
  }
}

template <class vertex>
graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+strlen(config)+1];
  char adjFile[strlen(iFile)+strlen(adj)+1];
  char idxFile[strlen(iFile)+strlen(idx)+1];
  *configFile = *adjFile = *idxFile = '\0';
  strcat(configFile,iFile);
  strcat(adjFile,iFile);
  strcat(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  long n;
  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
#ifdef WEIGHTED
  long m = size/(2*sizeof(uint));
#else
  long m = size/sizeof(uint);
#endif
  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();
  uintE* edges = (uintE*) s;

  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(intT)) { cout << "File size wrong\n"; abort(); }

  char* t = (char *) malloc(size);
  in3.read(t,size);
  in3.close();
  uintT* offsets = (uintT*) t;

  vertex* v = newA(vertex,n);
#ifdef WEIGHTED
  intE* edgesAndWeights = newA(intE,2*m);
  {parallel_for(long i=0;i<m;i++) {
    edgesAndWeights[2*i] = edges[i];
    edgesAndWeights[2*i+1] = edges[i+m];
    }}
  //free(edges);
#endif
  {parallel_for(long i=0;i<n;i++) {
    uintT o = offsets[i];
    uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors((uintE*)edges+o);
#else
      v[i].setOutNeighbors(edgesAndWeights+2*o);
#endif
  }}

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m);
#else
    intTriple* temp = newA(intTriple,m);
#endif
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);
#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif
    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);
    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);
    {parallel_for(long i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors((uintE*)inEdges+o);
#else
      v[i].setInNeighbors((intE*)(inEdges+2*o));
#endif
      }}
    free(tOffsets);
#ifndef WEIGHTED
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges);
    return graph<vertex>(v,n,m,mem);
#else
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights,inEdges);
    return graph<vertex>(v,n,m,mem);
#endif
  }
  free(offsets);
#ifndef WEIGHTED
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges);
  return graph<vertex>(v,n,m,mem);
#else
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights);
  return graph<vertex>(v,n,m,mem);
#endif
}

template <class vertex>
graph<vertex> readGraph(char* iFile, bool compressed, bool symmetric, bool binary, bool mmap) {
  if(binary) return readGraphFromBinary<vertex>(iFile,symmetric);
  else return readGraphFromFile<vertex>(iFile,symmetric,mmap);
}

template <class vertex>
graph<vertex> readCompressedGraph(char* fname, bool isSymmetric, bool mmap) {
  char* s;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    // Cannot mutate graph unless we copy.
    char *bytes = newA(char, S.n);
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    s = bytes;
  } else {
    ifstream in(fname,ifstream::in |ios::binary);
    in.seekg(0,ios::end);
    long size = in.tellg();
    in.seekg(0);
    cout << "size = " << size << endl;
    s = (char*) malloc(size);
    in.read(s,size);
    in.close();
  }

  long* sizes = (long*) s;
  long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  cout << "n = "<<n<<" m = "<<m<<" totalSpace = "<<totalSpace<<endl;
  cout << "reading file..."<<endl;

  uintT* offsets = (uintT*) (s+3*sizeof(long));
  long skip = 3*sizeof(long) + (n+1)*sizeof(intT);
  uintE* Degrees = (uintE*) (s+skip);
  skip+= n*sizeof(intE);
  uchar* edges = (uchar*)(s+skip);

  uintT* inOffsets;
  uchar* inEdges;
  uintE* inDegrees;
  if(!isSymmetric){
    skip += totalSpace;
    uchar* inData = (uchar*)(s + skip);
    sizes = (long*) inData;
    long inTotalSpace = sizes[0];
    cout << "inTotalSpace = "<<inTotalSpace<<endl;
    skip += sizeof(long);
    inOffsets = (uintT*) (s + skip);
    skip += (n+1)*sizeof(uintT);
    inDegrees = (uintE*)(s+skip);
    skip += n*sizeof(uintE);
    inEdges = (uchar*)(s + skip);
  } else {
    inOffsets = offsets;
    inEdges = edges;
    inDegrees = Degrees;
  }


  vertex *V = newA(vertex,n);
  parallel_for(long i=0;i<n;i++) {
    long o = offsets[i];
    uintT d = Degrees[i];
    V[i].setOutDegree(d);
    V[i].setOutNeighbors(edges+o);
  }

  if(sizeof(vertex) == sizeof(compressedAsymmetricVertex)){
    parallel_for(long i=0;i<n;i++) {
      long o = inOffsets[i];
      uintT d = inDegrees[i];
      V[i].setInDegree(d);
      V[i].setInNeighbors(inEdges+o);
    }
  }

  cout << "creating graph..."<<endl;
  Compressed_Mem<vertex>* mem = new Compressed_Mem<vertex>(V, s);

  graph<vertex> G(V,n,m,mem);
  return G;
}

/* 
    Populate a bitmap of where the hub vertices in a graph exist. 

    Also prepare maps from vertexID to index into local buffers
*/
template <class vertex>
void identifyHubs(graph<vertex> GA, Bitmap& bm, uintE*& idxMap, 
                  uintE*& inv_idxMap, size_t elemSz, uintE &hubCount,
                  bool isPageRank = false, bool useOutdeg = false)
{
    Timer t; 
    t.Start();
    auto numVertices = GA.n;
    vertex* Vertices = GA.V;
    pvector<uintT> degrees(numVertices);
    
    #pragma omp parallel for
    for (uintE v = 0; v < numVertices; ++v)
    {
        vertex vtx = Vertices[v];
        if (useOutdeg) 
            degrees[v] = vtx.getOutDegree();
        else
            degrees[v] = vtx.getInDegree();
    }
        
    __gnu_parallel::sort(degrees.begin(), degrees.end(), 
                         std::greater<uintT>());
    
    intT llcCap   = 1 << 25; //32MB conservative so that all duplicates of hubs fit in our 35MB LLC
    size_t addCost = (sizeof(uintE) + sizeof(bool));  //bitmap cost should be 1/8B per hub (best-case)
    intT numHubs {0};
    if (isPageRank)
        numHubs  = llcCap / ((omp_get_max_threads() * elemSz) + addCost);
    else
        numHubs  = llcCap / ((omp_get_max_threads() * elemSz) + sizeof(bool) + addCost);
    
    //assert(degrees[0] > degrees[numHubs] && degrees[numHubs] > degrees[numVertices-1]); //simple sanity check

    auto degThresh = degrees[numHubs]; //Hubs need to have degree greater than this

    /* Forming a map between original IDs and index into local buffer*/
    SlidingQueue<uint> queue(numHubs); 
    #pragma omp parallel
    {
        QueueBuffer<uintE> lqueue(queue); 
        #pragma omp for schedule(dynamic, 1024)
        for (uintE v = 0; v < numVertices; ++v)
        {
            if (useOutdeg)
            {
                if (Vertices[v].getOutDegree() > degThresh)
                {
                    bm.set_bit(v);
                    lqueue.push_back(v); 
                }
            }
            else
            {
                if (Vertices[v].getInDegree() > degThresh)
                {
                    bm.set_bit(v);
                    lqueue.push_back(v);
                }
            }
        }
        lqueue.flush();
    }
    queue.slide_window();
    
    hubCount   = queue.size();
    idxMap     = new uintE [GA.n]();
    inv_idxMap = new uintE [hubCount]();

    #pragma omp parallel for 
    for (auto q_iter = queue.begin(); q_iter < queue.end(); ++q_iter)
    {
        uintE vtxID = *q_iter;
        uintE posID = q_iter - queue.begin();
        idxMap[vtxID] = posID;
        inv_idxMap[posID] = vtxID;
    }
    t.Stop();
    t.PrintTime("HubDetect Time", t.Seconds());
}


/* 
  Preprocess a graph based on outdegrees or indegrees

  NOTE: This is NOT hub-sort. It is a proper sort based on degrees

  NOTE2: This is a relaxed version of degree-sort (the neighbor list of
         each vertex is NOT sorted : Different from `main-hubdup` results)
*/
template <class vertex>
graph<vertex> preprocessGraph(graph<vertex> GA, bool isSym, bool useOutdeg, 
                              pvector<uintE>& new_ids)
{
    Timer t; 
    t.Start();
    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    typedef std::pair<uintT, uintE> degree_nodeid_t; 
    pvector<degree_nodeid_t> degree_id_pairs(numVertices);
    if (!isSym) {
        /* directed graph */
        /* STEP I - collect degrees of all vertices */
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            vertex vtx = origG[v];
            if (useOutdeg) {
                degree_id_pairs[v] = std::make_pair(vtx.getOutDegree(), v);
            }
            else {
                degree_id_pairs[v] = std::make_pair(vtx.getInDegree(), v);
            }
        }

        /* Step II - sort the degrees in parallel */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(), 
                             std::greater<degree_nodeid_t>());

        /* Step III - make a remap based on the sorted degree list */
        pvector<uintT> degrees(numVertices);
        pvector<uintT> inv_degrees(numVertices);
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            degrees[v] = degree_id_pairs[v].first;
            auto origID = degree_id_pairs[v].second;
            new_ids[origID] = v;
            vertex vtx = origG[origID];
            if (useOutdeg) {
                inv_degrees[v] = vtx.getInDegree();
            }
            else {
                inv_degrees[v] = vtx.getOutDegree();
            }
        }
        //clearing space from degree pairs
        pvector<degree_nodeid_t>().swap(degree_id_pairs);

        /* Step IV - make a new vertex list for the new graph */
		pvector<uintT> offsets     = ParallelPrefixSum(degrees);
		pvector<uintT> inv_offsets = ParallelPrefixSum(inv_degrees);
        //clearing space from degree lists
        pvector<uintT>().swap(degrees);
        pvector<uintT>().swap(inv_degrees);
        #ifndef WEIGHTED
            uintE* outEdges = newA(uintE, numEdges);
            uintE* inEdges  = newA(uintE, numEdges);
        #else
            intE* outEdges = newA(intE, 2 * numEdges);
            intE* inEdges  = newA(intE, 2 * numEdges);
        #endif
        vertex* newV    = newA(vertex, numVertices);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (uintE v = 0; v < numVertices; ++v) {
            /* note that vertex IDs u and v belong to the space of original vertex IDs */
            //copy out-neighbors
            auto newID = new_ids[v];
            newV[newID].setOutDegree(origG[v].getOutDegree());
            #ifndef WEIGHTED
                if (useOutdeg)
                    newV[newID].setOutNeighbors(outEdges + offsets[newID]); 
                else
                    newV[newID].setOutNeighbors(outEdges + inv_offsets[newID]); 
            #else
                if (useOutdeg)
                    newV[newID].setOutNeighbors(outEdges + 2 * offsets[newID]); 
                else
                    newV[newID].setOutNeighbors(outEdges + 2 * inv_offsets[newID]); 
            #endif
            for (uintE u = 0; u < origG[v].getOutDegree(); ++u) {
                auto origNgh = origG[v].getOutNeighbor(u);
                newV[newID].setOutNeighbor(u, new_ids[origNgh]);
                #ifdef WEIGHTED
                    newV[newID].setOutWeight(u, origG[v].getOutWeight(u));
                #endif
            }
            
            //copy in-neighbors
            newV[newID].setInDegree(origG[v].getInDegree());
            #ifndef WEIGHTED
                if (useOutdeg) 
                    newV[newID].setInNeighbors(inEdges + inv_offsets[newID]); 
                else
                    newV[newID].setInNeighbors(inEdges + offsets[newID]); 
            #else
                if (useOutdeg) 
                    newV[newID].setInNeighbors(inEdges + 2 * inv_offsets[newID]); 
                else
                    newV[newID].setInNeighbors(inEdges + 2 * offsets[newID]); 
            #endif
            for (uintE u = 0; u < origG[v].getInDegree(); ++u) {
                auto origNgh = origG[v].getInNeighbor(u);
                newV[newID].setInNeighbor(u, new_ids[origNgh]);
                #ifdef WEIGHTED
                    newV[newID].setInWeight(u, origG[v].getInWeight(u));
                #endif
            }
        }

        /* Step V - make the new graph */ 
        Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,numVertices,numEdges,outEdges,inEdges);
        t.Stop();
        t.PrintTime("DegSort Time", t.Seconds());
        return graph<vertex>(newV,numVertices,numEdges,mem);
    }
    else {
        /* undirected graph */
        /* STEP I - collect degrees of all vertices */
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            vertex vtx = origG[v];
            degree_id_pairs[v] = std::make_pair(vtx.getOutDegree(), v);
        }

        /* Step II - sort the degrees in parallel */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(), 
                             std::greater<degree_nodeid_t>());

        /* Step III - make a remap based on the sorted degree list */
        pvector<uintT> degrees(numVertices);
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            degrees[v] = degree_id_pairs[v].first;
            auto origID = degree_id_pairs[v].second;
            new_ids[origID] = v;
        }
        //clearing space from degree pairs
        pvector<degree_nodeid_t>().swap(degree_id_pairs);

        /* Step IV - make a new vertex list for the new graph */
		pvector<uintT> offsets     = ParallelPrefixSum(degrees);
        //clearing space from degrees
        pvector<uintT>().swap(degrees);
        #ifndef WEIGHTED
            uintE* outEdges = newA(uintE, numEdges);
        #else
            intE* outEdges = newA(intE, 2 * numEdges);
        #endif
        vertex* newV    = newA(vertex, numVertices);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (uintE v = 0; v < numVertices; ++v) {
            /* note that vertex IDs u and v belong to the space of original vertex IDs */
            //copy neighbors
            auto newID = new_ids[v];
            newV[newID].setOutDegree(origG[v].getOutDegree());
            #ifndef WEIGHTED
                newV[newID].setOutNeighbors(outEdges + offsets[newID]); 
            #else
                newV[newID].setOutNeighbors(outEdges + 2 * offsets[newID]); 
            #endif
            for (uintE u = 0; u < origG[v].getOutDegree(); ++u) {
                auto origNgh = origG[v].getOutNeighbor(u);
                newV[newID].setOutNeighbor(u, new_ids[origNgh]);
                #ifdef WEIGHTED
                    newV[newID].setOutWeight(u, origG[v].getOutWeight(u));
                #endif
            }
        }

        /* Step V - make the new graph */ 
        Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,numVertices,numEdges,outEdges);
        t.Stop();
        t.PrintTime("DegSort Time", t.Seconds());
        return graph<vertex>(newV,numVertices,numEdges,mem);
    }
}

