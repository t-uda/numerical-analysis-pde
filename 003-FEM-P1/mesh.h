
#ifndef MESH_H_20170114204300
#define MESH_H_20170114204300

#include <stdbool.h>

#include "vector.h"

typedef struct vertex Vertex;
typedef struct edge Edge;
typedef struct triangle Triangle;

typedef struct vertex {
	size_t id;
	Vector2D pos;
	bool is_internal;
} Vertex;

typedef struct edge {
	Vertex * p;
	Vertex * q;
} Edge;

typedef struct triangle {
	Edge * edges[3];
	Vertex * vertices[3];
} Triangle;

typedef struct mesh {
	size_t nbv; // Number of vertices
	size_t nbe; // Number of edges
	size_t nbt; // Number of triangles
	Vertex * vertices;
	Edge * edges;
	Triangle * triangles;
} Mesh;

double prod_L2_tau(Triangle, Vertex *, Vertex *);

double semi_prod_H1_tau(Triangle, Vertex *, Vertex *);

Mesh read_mesh(char []);
void free_mesh(Mesh);

#endif

