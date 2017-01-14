
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"

// 行列式 |q-p, r-p|
double determinant(Vector2D p, Vector2D q, Vector2D r) {
	double a = q.x - p.x, b = r.x - p.x, c = q.y - p.y, d = r.y - p.y;
	return a * d - b * c;
}

// 三角形 tau の面積（の絶対値）の倍
double area_parallelogram(Triangle tau) {
	return fabs(determinant(tau.vertices[0]->pos, tau.vertices[1]->pos, tau.vertices[2]->pos));
}

// 三角形要素 tau 上での基底関数同士の L2 内積を返す
double prod_L2_tau(Triangle tau, Vertex * pi_ptr, Vertex * pj_ptr) {
	double d = area_parallelogram(tau);
	return (pi_ptr == pj_ptr) ? (d / 12.0) : (d / 24.0);
}

// 三角形 tau 上において頂点 *pi_ptr の対辺の内向き法線を返す
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
	Vector2D qr = vec2D_sub(q, r);
	Vector2D n = {-qr.y, qr.x};
	if (d < 0.0) { n.x = -n.x; n.y = -n.y; }
	return n;
}

// 三角形要素 tau 上での基底関数同士の H1 セミ内積（grad をかけた L2 内積）を返す
double semi_prod_H1_tau(Triangle tau, Vertex * pi_ptr, Vertex * pj_ptr) {
	double area = area_parallelogram(tau);
	Vector2D grad_i = inward_normal_vector(tau, pi_ptr); // area で割ると i 番目の基底関数の勾配
	Vector2D grad_j = inward_normal_vector(tau, pj_ptr); // j についても同様
	return 0.5 * vec2D_prod(grad_i, grad_j) / area; // 掛けて面積分 i.e. L2 内積
}

Mesh read_mesh(char mesh_file_name[]) {
	FILE * mesh_file = fopen(mesh_file_name, "r");
	if (mesh_file == NULL) {
		fprintf(stderr, "Failed to open the mesh file. ABORT.");
		exit(EXIT_FAILURE);
	}
	size_t nbv = 0; // Number of vertices
	size_t nbe = 0; // Number of edges
	size_t nbt = 0; // Number of triangles
	Vertex * vertices = NULL;
	Edge * edges = NULL;
	Triangle * triangles = NULL;
	// 一行ごとにメッシュファイルを読み込む
	char buffer[1024]; // ループ中，各行の文字列が buffer に読み込まれる
	while (fgets(buffer, sizeof(buffer), mesh_file) != NULL) {
		// まず nbv, nbe, nbt を読み込んで，それぞれの数の分だけメモリを確保する
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
			double x, y; // Vertex (id, x, y, is_internal)
			size_t p, q; // Edge (p, q)
			size_t e0, e1, e2; // Triangle (p, q, r, e0, e1, e2)
			int is_internal;
			if (sscanf(buffer, "Vertex %zu %lf %lf %d", &id, &x, &y, &is_internal) == 4) {
				// buffer から読み取ったデータを id 番目の頂点に設定する
				vertices[id].id = id;
				vertices[id].pos.x = x;
				vertices[id].pos.y = y;
				vertices[id].is_internal = is_internal;
			} else if (sscanf(buffer, "Edge %zu %zu %zu", &id, &p, &q) == 3) {
				// buffer から読み取ったデータを id 番目の辺に設定する
				edges[id].p = &vertices[p];
				edges[id].q = &vertices[q];
			} else if (sscanf(buffer, "Triangle %zu %zu %zu %zu", &id, &e0, &e1, &e2) == 4) {
				// buffer から読み取ったデータを id 番目の三角形に設定する
				Triangle * tau = &triangles[id];
				tau->edges[0] = &edges[e0];
				tau->edges[1] = &edges[e1];
				tau->edges[2] = &edges[e2];
				Vertex * p, * q;
				tau->vertices[0] = p = edges[e0].p;
				tau->vertices[1] = q = edges[e0].q;
				tau->vertices[2] = (edges[e1].p == p || edges[e1].p == q) ? edges[e1].q : edges[e1].p;
			}
		}
	}
	fclose(mesh_file); // 使い終わったメッシュファイルを閉じる
	Mesh mesh = {nbv, nbe, nbt, vertices, edges, triangles};
	return mesh;
}

void free_mesh(Mesh mesh) {
	free(mesh.triangles);
	free(mesh.edges);
	free(mesh.vertices);
}

