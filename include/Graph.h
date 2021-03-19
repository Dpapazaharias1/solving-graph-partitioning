/**@file  Graph.h
 *
 * @brief  Graph represesentation class consisting of the adjacency lists of its vertices
 *
 * @details This class reads either a standard edge type file or an adjacency lists
 * file (@see Graph::Graph for details). 
 *
 * @author Jose L. Walteros
 *
 * @version   0.0
 * @date      Feb 2020
 *
 *---+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/
#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <string>
#include <sstream>
#include <vector>

using namespace std;

class Graph
{

public:
      const char *name;       /**< Graph's name (filename) */
      int n;                  /**< Number of vertices in G */
      int m;                  /**< Number of edges in G*/
      int r;
      std::string graphtype;
      vector<int> EdgeTo;     /**< Vector that contains the adjacency list of each
                         * Vertex. The lists are appended one after the other */
      vector<int> EdgesBegin; /**< The position of the first neighbor of each vertex */
      vector<int> Degree;     /**< The degree of each vertex */
      vector<int> Weight;     /**<  The weight of each vertex */
      /**
 * Graph constructor: Receives a filename with the graph information to be read.
 * The expected format is either the edge list or the adjacency lists of the graph.
 *
 * The vertices are expected to be
 * labeled from 0 to n-1.
 *
 * @param[in] filename : Name of the file with the graph information.
 */
      Graph(
          const char *filetype, // -e for edge list -a for adjacency list.
          const char *filename, // Input file either: 1. adjacency list, or 2. edge list
          const char *weighted  // Input file is weighted or not
      );

      /**
 * SP: Calculates the shortest path from a source sacording to a given cost vector
 *
 * The costs are expected to be a vector similar to EdgeTo, but having the cost in each position instead of the vertex
 *
 * @param[in] CostTo : const information.
 * @param[in] s : source.
 * @param[in] parent : parent of each node.
 * @param[in] cost : cost of SP to each node from s.
 */
      void SP(
          std::vector<double> CostTo,
          int s,
          double *parent,
          double *cost);

      /**
 * Prim: Find a minimum cost cover tree of size r + 1 starting from root s 
 * Parent list will contain the edge its connected to
 * @param[in] CostTo: const information
 * @param[in] s : source
 * @param[in] r : maximum component weight
 * @param[in] parent : parent of each node.
 * @param[in] cost : cost of edge in MST
 * 
 **/

      std::vector<int> Prim(
          std::vector<double> CostTo,
          int s,
          int r,
          double *parent,
          double *cost);

      /**
 * Prints the full description of the graph:
 *
 * <name n m delta Delta>
 * <v(name): adjacency list(v)>
 *
 * For example:
 *
 * test.graph 10 13 1 5
 * 0(4): 1 2 3 9
 * 1(2): 0 2
 * 2(2): 0 1
 * 3(5): 0 4 5 7 8
 * 4(3): 3 5 6
 * 5(4): 3 4 6 8
 * 6(2): 4 5
 * 7(1): 3
 * 8(2): 3 5
 * 9(1): 0
 */
      void print();
      /**
 * Prints a short description of graph:
 *
 * <name n m delta Delta>
 *
 * For example:
 *
 * test.graph 10 13 1 5
 */
      void printShort();
};

#endif
