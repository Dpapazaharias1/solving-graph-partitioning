/**@file  Graph.cpp
 *
 * @brief  Graph represesentation class consisting of the adjacency lists of its vertices
 *
 * @details This class reads either a standard edge type file or an adjacency lists
 * file (@see Graph::Graph for details). The Graph class also implements a procedure
 * that generates a degeneracy ordering of the given graph (@see Graph::degeneracyOrdering)
 * and a procedure that outputs the complement of the subgraph induced by a given
 * vertex v and it's neighbors to the right in the degeneracy ordering.
 *
 * @author Jose L. Walteros
 *
 * @version   0.0
 * @date      March 2017
 *
 *---+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <set>
#include <limits.h>
#include <stdlib.h>
#include <algorithm>
#include "Graph.h"
#include "MinHeap.h"

Graph::Graph(
	const char *filetype,
	const char *filename,
	const char *weighted)
{
	name = filename;

	ifstream file(filename);
	if (!file)
	{
		cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
	}

	file >> n >> m >> r >> graphtype;

	if (n * m == 0)
	{
		cerr << "ERROR: when reading the graph from file '" << filename;
	}

	else
	{
		/**
                 * The input files have several inconsistencies including duplicated edges.
                 * The adjacency lists are first filtered by a set that removes duplicates.
                 * Additionally, this filtering process sorts the adjacency lists, which is
                 * quite useful. The filtering process increases the time required to generate
                 * the graph, but removes duplicates. Cannot be removed unless the input files
                 * are fixed first.
                 */
		int i, j, w, counter;
		int real_m = 0;
		EdgesBegin = vector<int>(n, 0);
		Degree = vector<int>(n, 0);
		Weight = vector<int>(n, 1);
		vector<set<int>> adjLists(n, set<int>());
		/**
                 * If nodes are weighted we randomly select node weights
                 * In addition we update r to reflect the proportion of (r/n)->(r_w / weight(V))
                 **/

		std::cout << "n, r: " << n << ", " << r << std::endl;

		total_weight = n;
		bool is_weighted = (strcmp(weighted, "-w") == 0);
		if (is_weighted)
		{
			srand(1);
			for (int v = 0; v < n; v++)
			{
				w = rand() % 1000 + 1;
				Weight[v] = w;
				total_weight += w;
			}
			float r_pct = (float)r / (float)n;
			std::cout << r_pct << std::endl;
			r = (int)(total_weight * r_pct);
		}

		std::cout << "w(G), r: " << total_weight << ", " << r << std::endl;
		//std::cout << "weight" << std::endl;
		/**
                 * Procedure for reading edge list input files. The adjacency lists are first
                 * created as sets to ensure there are no duplicated edges and loops.
                 */
		if (strcmp(filetype, "-e") == 0)
		{
			for (int e = 0; e < m; e++)
			{
				file >> i >> j;
				if (i != j)
				{
					adjLists[i].insert(j).second;
					adjLists[j].insert(i).second;
					Degree[i]++;
					Degree[j]++;
					real_m++;
				}
			}
		}
		else
		{
			if (strcmp(filetype, "-a") == 0) // Procedure for reading edge list input files
			{
				string str;
				i = 0;
				std::getline(file, str);
				while (std::getline(file, str))
				{
					istringstream ss(str);
					while (ss >> j)
					{
						if (i != j)
						{
							adjLists[i].insert(j).second;
							Degree[i]++;
						}
					}
					real_m += Degree[i];
					i++;
				}
				real_m = real_m / 2;
			}
		}
		m = real_m;
		EdgeTo = vector<int>(2 * m);
		EdgeCost = vector<int>(2 * m, 1);
		counter = 0;
		/**
                 * The following procedure extracts the adjacency lists from the filter created by the sets. We dont keep the sets
                 * because vectors tend to perform better. The extra I/O overhead of the filter is a result of the inconsistencies in most
                 * of the graph files used for testing.
                 */
		srand(1000);
		for (int i = 0; i < n; i++)
		{
			EdgesBegin[i] = counter;

			for (std::set<int>::iterator it = adjLists[i].begin(); it != adjLists[i].end(); ++it)
			{
				if (is_weighted)
				{
					EdgeCost[counter] = rand() % 1000 + 1;
				}
				EdgeTo[counter++] = *it;
			}
		}

		if (is_weighted)
		{
			for (int i = 0; i < n; i++)
			{
				for (int k = EdgesBegin[i]; k < EdgesBegin[i] + Degree[i]; k++)
				{
					int j = EdgeTo[k];
					if (i < j)
					{
						for (int l = EdgesBegin[j]; l < EdgesBegin[j] + Degree[j]; l++)
						{
							if (i == EdgeTo[l])
							{
								EdgeCost[l] = EdgeCost[k];
							}
						}
					}
				}
			}
		}
	}
}

void Graph::SP(
	std::vector<double> CostTo,
	int s,
	int *parent,
	double *cost)
{
	for (int i = 0; i < n; i++)
	{
		parent[i] = -1;
		cost[i] = INT_MAX;
	}
	MinHeap Q(n);
	for (int i = 0; i < n; i++)
	{
		Q.array[i] = MinHeapNode(i, cost[i]);
		Q.pos[i] = i;
		Q.size++;
	}
	cost[s] = 0;
	Q.pos[s] = 0;
	Q.pos[0] = s;
	Q.array[s].key = 0;
	MinHeapNode temp = Q.array[s];
	Q.array[s] = Q.array[0];
	Q.array[0] = temp;
	while (!Q.isEmpty())
	{
		int min = Q.extractMin().v;
		for (int j = EdgesBegin[min]; j < Degree[min] + EdgesBegin[min]; j++)
		{
			double temp = cost[min] + CostTo[j];
			if (temp < cost[EdgeTo[j]])
			{
				cost[EdgeTo[j]] = temp;
				parent[EdgeTo[j]] = min;
				Q.decreaseKey(EdgeTo[j], temp);
			}
		}
	}
}

std::vector<int> Graph::Prim(
	std::vector<double> CostTo,
	int s,
	int r,
	int *parent,
	double *cost)
{
	for (int i = 0; i < n; i++)
	{
		parent[i] = -1;
		cost[i] = INT_MAX;
	}
	MinHeap Q(n);
	for (int i = 0; i < n; i++)
	{
		Q.array[i] = MinHeapNode(i, cost[i]);
		Q.pos[i] = i;
		Q.size++;
	}
	cost[s] = 0;
	Q.pos[s] = 0;
	Q.pos[0] = s;
	Q.array[s].key = 0;
	MinHeapNode temp = Q.array[s];
	Q.array[s] = Q.array[0];
	Q.array[0] = temp;

	std::vector<int> treeNodes;
	int weight = 0;

	//std::vector<bool> is_leaf(n, false);
	std::vector<bool> in_tree(n, false);
	in_tree[s] = true;
	//is_leaf[s] = true;
	//std::vector<int> tree_deg(n, 0);

	while (!Q.isEmpty() && weight < r + 1)
	{
		int min = Q.extractMin().v;
		in_tree[min] = true;
		treeNodes.push_back(min);
		weight += Weight[min]; // TODO: Need to change to vertex weight
		for (int j = EdgesBegin[min]; j < Degree[min] + EdgesBegin[min]; j++)
		{
			double temp = CostTo[j];
			if (temp < cost[EdgeTo[j]] && !in_tree[EdgeTo[j]])
			{
				cost[EdgeTo[j]] = temp;
				for(int k = EdgesBegin[EdgeTo[j]]; k <  EdgesBegin[EdgeTo[j]] + Degree[EdgeTo[j]]; k++)
				{
					if (EdgeTo[k] == min)
					{
						parent[EdgeTo[j]] = k;
					}
				}
				//parent[EdgeTo[j]] = j;
				Q.decreaseKey(EdgeTo[j], temp);
			}
			
		}
		/*
		tree_deg[min]++;
		is_leaf[min] = true;
		if (parent[min] > -1)
		{
			int j = EdgeTo[parent[min]];
			tree_deg[j]++;
			if (is_leaf[j])
			{
				is_leaf[j] = false;
			}
		}
		*/
		
		
	}
	/*
	std::vector<int> leaf_nodes;
	int min_weight = r;
	int min_leaf = 0;
	for(int i = 0; i < n; i++)
	{	
		if (is_leaf[i] || tree_deg[i] == 1)
		{
			std::cout << "Leaf " << i << ", weight:" << Weight[i] << std::endl;
			leaf_nodes.push_back(i);
			if (Weight[i] < min_weight)
			{
				min_weight = Weight[i];
				min_leaf = leaf_nodes.size() - 1;
			}
		}
	}
	
	//std::cout << "Leaf with min weight " << leaf_nodes[min_leaf] << ", " <<  min_weight << std::endl; 
	while( weight - min_weight > r)
	{
		
		int leaf = leaf_nodes[min_leaf];
		int p = EdgeTo[parent[leaf]];
		tree_deg[p]--;
		in_tree[leaf] = false;
		leaf_nodes.erase(leaf_nodes.begin() + min_leaf);
		parent[leaf_nodes[min_leaf]] = -1;
		weight -= min_weight;
		min_weight = r;
		if(tree_deg[p] == 1)
		{
			is_leaf[p] = true;
			leaf_nodes.push_back(p);
		}
		for( int i = 0; i < leaf_nodes.size(); i++)
		{
			int j = leaf_nodes[i];
			if (Weight[j] < min_weight)
			{
				min_weight = Weight[j];
				min_leaf = i;
			}
		}
		//std::cout << "Tree not a minimal cover, remove node " << leaf << std::endl;
		//std::cout << "New weight: " << weight << ", new min leaf: " << min_weight<<std::endl;
	}
	for(int i = 0; i < n; i++)
	{
		if( in_tree[i])
		{
			treeNodes.push_back(i);
		}
	}
	*/
	return treeNodes; 
}

void Graph::print()
{
	printShort();
	for (int i = 0; i < n; i++)
	{
		std::cout << i << "(" << Degree[i] << "): ";
		for (int j = EdgesBegin[i]; j < Degree[i] + EdgesBegin[i]; j++)
			std::cout << EdgeTo[j] << " ";
		std::cout << "\n";
	}
	std::cout << "\n";
}

void Graph::printShort()
{
	cout << " " << name << " " << n << " " << m << "\n";
}
