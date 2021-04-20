#ifndef _FORMULATION_H_
#define _FORMULATION_H_

void triangle (const Graph &G);
void flow(Graph &G, int p);
void path(Graph &G, int p);

std::vector<double> PathFormulation(Graph &G);

#endif