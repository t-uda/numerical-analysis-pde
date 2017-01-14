// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct vector2d Vector2D;
typedef struct vertex Vertex;
typedef struct edge Edge;
typedef struct triangle Triangle;

typedef struct vector2d {
	double x;
	double y;
} Vector2D;

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

Vector2D vec_sub(Vector2D a, Vector2D b) {
	Vector2D c = {a.x - b.x, a.y - b.y};
	return c;
}

double vec_prod(Vector2D a, Vector2D b) {
	return a.x * b.x + a.y * b.y;
}

double length(Vector2D v) {
	return hypot(v.x, v.y);
}

// // // // // // // // // * キ ケ ン * // // // // // // // // // //
// // // // // malloc 関数の使用には十分注意すること！ // // // // //
// メモリ確保命令 malloc, calloc を使い浮動小数点数型の配列を動的に
// 確保する．使い終わったら必ず free() すること．
void check_allocated_or_abort(void * ptr) {
	if (ptr == NULL) {
		fprintf(stderr, "Failed to allocate a vector. ABORT.\n");
		exit(EXIT_FAILURE);
	}
}
double * allocate_real_vector(size_t n) {
	double * vector = calloc(n, sizeof(double));
	check_allocated_or_abort(vector);
	return vector;
}
// // // // // // // // // * キ ケ ン * // // // // // // // // // //

// 行列式 |q-p, r-p|
double determinant(Vector2D p, Vector2D q, Vector2D r) {
	double a = q.x - p.x, b = r.x - p.x, c = q.y - p.y, d = r.y - p.y;
	return a * d - b * c;
}

double area_parallelogram(Triangle tau) {
	return fabs(determinant(tau.vertices[0]->pos, tau.vertices[1]->pos, tau.vertices[2]->pos));
}

// 三角形要素 tau 上での基底関数同士の L2 内積を返す
double prod_L2_tau(Triangle tau, size_t i, size_t j) {
	double d = area_parallelogram(tau);
	return (i == j) ? (d / 12.0) : (d / 24.0);
}

// i 番目の頂点の対辺の内向き法線を返す
Vector2D inward_normal_vector(Triangle tau, Vertex * pi_ptr) {
	Vector2D p, q, r;
	p = pi_ptr->pos;
	if (tau.vertices[0] == pi_ptr) {
		q = tau.vertices[1]->pos;
		r = tau.vertices[2]->pos;
	} else if (tau.vertices[1] == pi_ptr) {
		q = tau.vertices[2]->pos;
		r = tau.vertices[0]->pos;
	} else if (tau.vertices[2] == pi_ptr) {
		q = tau.vertices[0]->pos;
		r = tau.vertices[1]->pos;
	}
	double d = determinant(p, q, r);
	Vector2D qr = vec_sub(q, r);
	Vector2D n = {-qr.y, qr.x};
	if (d < 0.0) { n.x = -n.x; n.y = -n.y; }
	return n;
}

// 三角形要素 tau 上での H1 セミ内積（grad をかけた L2 内積）を返す
double semi_prod_H1_tau(Triangle tau, Vertex * pi_ptr, Vertex * pj_ptr) {
	double area = area_parallelogram(tau);
	Vector2D grad_i = inward_normal_vector(tau, pi_ptr); // area で割ると i 番目の基底関数の勾配
	Vector2D grad_j = inward_normal_vector(tau, pj_ptr); // j についても同様
	return 0.5 * vec_prod(grad_i, grad_j) / area; // 掛けて面積分 i.e. L2 内積
}

int main(int argc, char * argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Specify a mesh file.");
		return EXIT_FAILURE;
	}
	FILE * mesh_file = fopen(argv[1], "r");
	if (mesh_file == NULL) {
		fprintf(stderr, "Failed to open the mesh file.");
		return EXIT_FAILURE;
	}
	size_t nbv = 0;
	size_t nbe = 0;
	size_t nbt = 0;
	Vertex * vertices = NULL;
	Edge * edges = NULL;
	Triangle * triangles = NULL;
	char buffer[1024];
	while (fgets(buffer, sizeof(buffer), mesh_file) != NULL) {
		if (nbv == 0 || nbe == 0 || nbt == 0) {
			if (sscanf(buffer, "NbVertices %zu", &nbv) == 1) {
				vertices = malloc(nbv * sizeof(Vertex));
				check_allocated_or_abort(vertices);
			} else if (sscanf(buffer, "NbEdges %zu", &nbe) == 1) {
				edges = malloc(nbe * sizeof(Edge));
				check_allocated_or_abort(edges);
			} else if (sscanf(buffer, "NbTriangles %zu", &nbt) == 1) {
				triangles = malloc(nbt * sizeof(Triangle));
				check_allocated_or_abort(triangles);
			}
		} else {
			size_t id;
			double x, y; // Vertex (x, y)
			size_t p, q; // Edge (p, q)
			size_t e0, e1, e2; // Triangle (e1, e2, e3)
			int is_internal;
			if (sscanf(buffer, "Vertex %zu %lf %lf %d", &id, &x, &y, &is_internal) == 4) {
				vertices[id].id = id;
				vertices[id].pos.x = x;
				vertices[id].pos.y = y;
				vertices[id].is_internal = is_internal;
			//	printf("Vertex#%zu (%f, %f)\n", id, x, y);
			} else if (sscanf(buffer, "Edge %zu %zu %zu", &id, &p, &q) == 3) {
				edges[id].p = &vertices[p];
				edges[id].q = &vertices[q];
			//	printf("Edge#%zu (%zu, %zu)\n", id, p, q);
			} else if (sscanf(buffer, "Triangle %zu %zu %zu %zu", &id, &e0, &e1, &e2) == 4) {
				Triangle * tau = &triangles[id];
				tau->edges[0] = &edges[e0];
				tau->edges[1] = &edges[e1];
				tau->edges[2] = &edges[e2];
				Vertex * p, * q;
				tau->vertices[0] = p = edges[e0].p;
				tau->vertices[1] = q = edges[e0].q;
				tau->vertices[2] = (edges[e1].p == p || edges[e1].p == q) ? edges[e1].q : edges[e1].p;
			//	printf("Triangle#%zu (%zu, %zu, %zu)\n", id, e0, e1, e2);
			}
		}
	}
	fclose(mesh_file);

	size_t nbu = 0;
	size_t * internal_ids = malloc(nbv * sizeof(size_t));
	check_allocated_or_abort(internal_ids);
	Vertex ** internal_vertices = malloc(nbv * sizeof(Vertex *));
	check_allocated_or_abort(internal_vertices);
	for (size_t m = 0; m < nbv; m++) {
		if (vertices[m].is_internal) {
			size_t id = nbu++;
			internal_ids[m] = id;
			internal_vertices[id] = &vertices[m];
		}
	}
	double * a_data = allocate_real_vector(nbu * nbu);
	double ** a = malloc(nbu * sizeof(double **));
	check_allocated_or_abort(a);
	for (size_t i = 0; i < nbu; i++) {
		a[i] = &a_data[i * nbu];
	}
	for (size_t e_id = 0; e_id < nbe; e_id++) {
		Edge * e = &edges[e_id];
		if (!e->p->is_internal || !e->q->is_internal) continue;
		size_t i = internal_ids[e->p->id], j = internal_ids[e->q->id];
		for (size_t k = 0; k < nbt; k++) {
			Triangle tau = triangles[k];
			if (e == tau.edges[0] || e == tau.edges[1] || e == tau.edges[2]) {
				a[i][j] += semi_prod_H1_tau(tau, e->p, e->q);
			}
		}
		a[j][i] = a[i][j];
	}
	for (size_t i = 0; i < nbu; i++) {
		Vertex * p = internal_vertices[i];
		for (size_t k = 0; k < nbt; k++) {
			Triangle tau = triangles[k];
			if (p == tau.vertices[0] || p == tau.vertices[1] || p == tau.vertices[2]) {
				a[i][i] += semi_prod_H1_tau(tau, p, p);
			}
		}
	}
	for (size_t i = 0; i < nbu; i++) {
		for (size_t j = i; j < nbu; j++) {
			fprintf(stderr, "a[%zu][%zu] = %le\n", i, j, a[i][j]);
		}
	}
	free(a);
	free(a_data);
	free(internal_ids);
	free(triangles);
	free(edges);
	free(vertices);
	return 0;
}

