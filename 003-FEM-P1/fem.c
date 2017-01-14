// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

double f(double x, double y) {
//	return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
	return 32.0 * (x * (1.0 - x) + y * (1.0 - y));
}

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

// returns the sum a+b
Vector2D vec2D_sub(Vector2D a, Vector2D b) {
	Vector2D c = {a.x - b.x, a.y - b.y};
	return c;
}

// returns the inner product a.b
double vec2D_prod(Vector2D a, Vector2D b) {
	return a.x * b.x + a.y * b.y;
}

// returns the length |v|
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

void solve_CG(size_t, double *, double **, double *);

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

	size_t nbu = 0; // Number of unknowns （＝節点の数＝内部頂点の数） < nbv
	size_t * internal_ids = malloc(nbv * sizeof(size_t)); // 頂点番号を添え字とするその節点番号の配列
	check_allocated_or_abort(internal_ids);
	Vertex ** internal_vertices = malloc(nbv * sizeof(Vertex *)); // 節点番号を添え字とするその頂点へのポインタの配列
	check_allocated_or_abort(internal_vertices);
	for (size_t m = 0; m < nbv; m++) {
		if (vertices[m].is_internal) {
			size_t id = nbu++;
			internal_ids[m] = id;
			internal_vertices[id] = &vertices[m];
		}
	}

	// サイズ nbu*nbu の double 型配列を確保し，それを行列として使う
	double * a_data = allocate_real_vector(nbu * nbu);
	double ** a = malloc(nbu * sizeof(double **));
	check_allocated_or_abort(a);
	for (size_t i = 0; i < nbu; i++) {
		a[i] = &a_data[i * nbu];
	}

	// a[i][j] (i != j) を計算する．
	// i 番節点と j 番節点を結ぶ辺があるところが非零要素であるから，辺でループを回す
	for (size_t e_id = 0; e_id < nbe; e_id++) {
		Edge * e = &edges[e_id];
		if (!e->p->is_internal || !e->q->is_internal) continue; // 辺の端点が境界にある場合は考えなくてよい
		size_t i = internal_ids[e->p->id], j = internal_ids[e->q->id];
		// 辺 e をもつ三角形 tau の上で内積を足しあげる
		for (size_t k = 0; k < nbt; k++) {
			Triangle tau = triangles[k];
			if (e == tau.edges[0] || e == tau.edges[1] || e == tau.edges[2]) {
				a[i][j] += semi_prod_H1_tau(tau, e->p, e->q);
			}
		}
		// 対称性
		a[j][i] = a[i][j];
	}
	// 対角要素 a[i][i] を計算する
	for (size_t i = 0; i < nbu; i++) {
		Vertex * p = internal_vertices[i];
		for (size_t k = 0; k < nbt; k++) {
			Triangle tau = triangles[k];
			if (p == tau.vertices[0] || p == tau.vertices[1] || p == tau.vertices[2]) {
				a[i][i] += semi_prod_H1_tau(tau, p, p);
			}
		}
	}

	// 右辺ベクトル b を計算する
	double * b = allocate_real_vector(nbu);
	for (size_t i = 0; i < nbu; i++) {
		Vertex * p = internal_vertices[i];
		// 節点 p を頂点にもつ三角形 tau について積分を足しあげる
		for (size_t k = 0; k < nbt; k++) {
			Triangle tau = triangles[k];
			if (p == tau.vertices[0] || p == tau.vertices[1] || p == tau.vertices[2]) {
				// ソースターム f の P1 要素補間 f_h と基底関数の L2 内積を計算している：
				//     f_h = \sum_q f(q) \varphi_q,
				//     < f_h, \varphi_p > = \sum_q f(q) < \varphi_q, \varphi_p >.
				for (size_t v = 0; v < 3; v++) {
					Vertex * q = tau.vertices[v];
					b[i] += f(q->pos.x, q->pos.y) * prod_L2_tau(tau, p, q);
				}
			}
		}
	}

	// 方程式 a * x = b を CG 法で解く
	double * x = allocate_real_vector(nbu);
	solve_CG(nbu, x, a, b);

	// 計算結果を出力する
	for (size_t m = 0; m < nbv; m++) {
		Vertex p = vertices[m];
		double u = 0.0;
		if (p.is_internal) {
			size_t i = internal_ids[m];
			u = x[i];
		}
		printf("%1.17le\t%1.17le\t%1.17le\n", p.pos.x, p.pos.y, u);
	}

	free(x);
	free(b);
	free(a);
	free(a_data);
	free(internal_ids);
	free(triangles);
	free(edges);
	free(vertices);
	return 0;
}

// x := y
void vec_set(size_t n, double * x, double * y) {
	for (size_t i = 0; i < n; i++) {
		x[i] = y[i];
	}
}

// x += y
void vec_add(size_t n, double * x, double * y) {
	for (size_t i = 0; i < n; i++) {
		x[i] += y[i];
	}
}

// x -= y
void vec_sub(size_t n, double * x, double * y) {
	for (size_t i = 0; i < n; i++) {
		x[i] -= y[i];
	}
}

// x *= r componentwise
void vec_mul(size_t n, double r, double * x) {
	for (size_t i = 0; i < n; i++) {
		x[i] *= r;
	}
}

// res = A * x
void mat_vec_prod(size_t n, double * res, double ** a, double * x) {
	for (size_t i = 0; i < n; i++) {
		res[i] = 0.0;
		for (size_t j = 0; j < n; j++) {
			res[i] += a[i][j] * x[j];
		}
	}
}

// returns the inner product x.y
double vec_prod(size_t n, double * x, double * y) {
	double res = 0.0;
	for (size_t i = 0; i < n; i++) {
		res += x[i] * y[i];
	}
	return res;
}

// returns the bilinear form A[x, y] = x.(A*y)
double bilinear_form(size_t n, double * x, double ** a, double * y) {
	double * ay = allocate_real_vector(n);
	mat_vec_prod(n, ay, a, y);
	double xay = vec_prod(n, x, ay);
	free(ay);
	return xay;
}

// 線形方程式 A x = b を共役勾配法 (conjugate gradient method) で解く
// （ただし行列 A は正定値対称行列であることを仮定）
void solve_CG(size_t n, double * x, double ** a, double * b) {
	double * tmp = allocate_real_vector(n);
	double * r = allocate_real_vector(n);
	double * p = allocate_real_vector(n);
	double * alpha_p = allocate_real_vector(n);
	// r = b - A x
	mat_vec_prod(n, tmp, a, x);
	vec_add(n, r, b);
	vec_sub(n, r, tmp);
	// p = r
	vec_add(n, p, r);
	for (size_t iter = 0; iter < n; iter++) {
		// alpha = r.r / A[p, p]
		double r_sqr = vec_prod(n, r, r);
		double alpha = r_sqr / bilinear_form(n, p, a, p);
		// alpha_p = alpha * p
		vec_set(n, alpha_p, p);
		vec_mul(n, alpha, alpha_p);
		// x += alpha_p
		vec_add(n, x, alpha_p);
		// r -= A * alpha_p
		mat_vec_prod(n, tmp, a, alpha_p);
		vec_sub(n, r, tmp);
		if (vec_prod(n, r, r) < 1.0e-14) break;
		double beta = vec_prod(n, r, r) / r_sqr;
		// p := r + beta * p
		vec_set(n, tmp, p);
		vec_set(n, p, r);
		vec_mul(n, beta, tmp);
		vec_add(n, p, tmp);
	}
	free(alpha_p);
	free(p);
	free(r);
	free(tmp);
}

