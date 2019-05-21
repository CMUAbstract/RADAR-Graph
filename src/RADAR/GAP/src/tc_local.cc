// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <cassert>
#include <omp.h>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"


/*
GAP Benchmark Suite
Kernel: Triangle Counting (TC)
Author: Scott Beamer

Will count the number of triangles (cliques of size 3)

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

This implementation reduces the search space by counting each triangle only
once. A naive implementation will count the same triangle six times because
each of the three vertices (u, v, w) will count it in both ways. To count
a triangle only once, this implementation only counts a triangle if u > v > w.
Once the remaining unexamined neighbors identifiers get too big, it can break
out of the loop, but this requires that the neighbors to be sorted.

Another optimization this implementation has is to relabel the vertices by
degree. This is beneficial if the average degree is high enough and if the
degree distribution is sufficiently non-uniform. To decide whether or not
to relabel the graph, we use the heuristic in WorthRelabelling.
*/


using namespace std;

size_t OrderedCount(const Graph &g) {
  Timer tm;
  tm.Start();
  size_t total = 0;
  pvector<size_t> triangles(g.num_nodes(), 0);
  
  //create thread private storage
  NodeID llcCap   = 1 << 25; //32MB conservative so that all duplicates of hubs fit in our 35MB LLC
  NodeID numHubs  = llcCap / (omp_get_max_threads() * sizeof(size_t));
  size_t** triangles_dup = new size_t* [omp_get_max_threads()];
  assert(triangles_dup != nullptr);
  for (int tid = 0; tid < omp_get_max_threads(); ++tid)
  {
    triangles_dup[tid] = new size_t [numHubs]();
    assert(triangles_dup[tid] != nullptr);
  }

  #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    for (NodeID v : g.out_neigh(u)) {
      if (v > u)
        break;
      auto it = g.out_neigh(u).begin();
      for (NodeID w : g.out_neigh(v)) {
        if (w > v)
          break;
        while (*it < w)
          it++;
        if (w == *it)
        {
          total++;
          if (u < numHubs)
            ++triangles_dup[omp_get_thread_num()][u];
          else
            fetch_and_add(triangles[u], 1);

          if (v < numHubs)
            ++triangles_dup[omp_get_thread_num()][v];
          else
            fetch_and_add(triangles[v], 1);
        
          if (w < numHubs)
            ++triangles_dup[omp_get_thread_num()][w];
          else
            fetch_and_add(triangles[w], 1);
        }
      }
    }
  }

  //merge thread private duplicates
  #pragma omp parallel for schedule(static, 64 / sizeof(size_t))
  for (NodeID u = 0; u < numHubs; ++u)
  {
    for (int tid = 0; tid < omp_get_max_threads(); ++tid)
    {
      triangles[u] += triangles_dup[tid][u];
    }
  }
  
  //free memory
  for (int tid = 0; tid < omp_get_max_threads(); ++tid)
  {
    delete[] triangles_dup[tid];
  }
  delete[] triangles_dup;

  tm.Stop();
  std::cout << "Run Time (in s) = " << tm.Seconds() << std::endl;

  //verify output
  size_t tricnt {0};
  #pragma omp parallel for reduction(+ : tricnt)
  for (NodeID u = 0; u < g.num_nodes(); ++u)
  {
    tricnt += triangles[u];
  }
  assert(tricnt == 3 * total);
  return total;
}


// heuristic to see if sufficently dense power-law graph
bool WorthRelabelling(const Graph &g) {
  int64_t average_degree = g.num_edges() / g.num_nodes();
  if (average_degree < 10)
    return false;
  SourcePicker<Graph> sp(g);
  int64_t num_samples = min(int64_t(1000), g.num_nodes());
  int64_t sample_total = 0;
  pvector<int64_t> samples(num_samples);
  for (int64_t trial=0; trial < num_samples; trial++) {
    samples[trial] = g.out_degree(sp.PickNext());
    sample_total += samples[trial];
  }
  sort(samples.begin(), samples.end());
  double sample_average = static_cast<double>(sample_total) / num_samples;
  double sample_median = samples[num_samples/2];
  return sample_average / 1.3 > sample_median;
}


// uses heuristic to see if worth relabeling
size_t Hybrid(const Graph &g) {
  if (WorthRelabelling(g))
    std::cout << "Graph supposed to be reordered\n";
  else
    std::cout << "Graph NOT supposed to be reordered\n";

  return OrderedCount(Builder::RelabelByDegree(g));
}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
  cout << total_triangles << " triangles" << endl;
}


// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total) {
  size_t total = 0;
  vector<NodeID> intersection;
  intersection.reserve(g.num_nodes());
  for (NodeID u : g.vertices()) {
    for (NodeID v : g.out_neigh(u)) {
      auto new_end = set_intersection(g.out_neigh(u).begin(),
                                      g.out_neigh(u).end(),
                                      g.out_neigh(v).begin(),
                                      g.out_neigh(v).end(),
                                      intersection.begin());
      intersection.resize(new_end - intersection.begin());
      total += intersection.size();
    }
  }
  total = total / 6;  // each triangle was counted 6 times
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "triangle count");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  BenchmarkKernel(cli, g, Hybrid, PrintTriangleStats, TCVerifier);
  return 0;
}
