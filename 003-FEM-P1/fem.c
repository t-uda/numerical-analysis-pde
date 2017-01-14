// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"
#include "mesh.h"

double f(double x, double y) {
//	return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
//	return 32.0 * (x * (1.0 - x) + y * (1.0 - y));
	return 1.0;
}

int main(int argc, char * argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Specify a mesh file.");
		exit(EXIT_FAILURE);
	}
	Mesh mesh = read_mesh(argv[1]);

	size_t nbu = 0; // Number of unknowns （＝節点の数＝内部頂点の数） < mesh.nbv
	size_t * internal_ids = malloc(mesh.nbv * sizeof(size_t)); // 頂点番号を添え字とするその節点番号の配列
	check_allocated_or_abort(internal_ids);
	Vertex ** internal_vertices = malloc(mesh.nbv * sizeof(Vertex *)); // 節点番号を添え字とするその頂点へのポインタの配列
	check_allocated_or_abort(internal_vertices);
	for (size_t m = 0; m < mesh.nbv; m++) {
		if (mesh.vertices[m].is_internal) {
			size_t id = nbu++;
			internal_ids[m] = id;
			internal_vertices[id] = &mesh.vertices[m];
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
	for (size_t e_id = 0; e_id < mesh.nbe; e_id++) {
		Edge * e = &mesh.edges[e_id];
		if (!e->p->is_internal || !e->q->is_internal) continue; // 辺の端点が境界にある場合は考えなくてよい
		size_t i = internal_ids[e->p->id], j = internal_ids[e->q->id];
		// 辺 e をもつ三角形 tau の上で内積を足しあげる
		for (size_t k = 0; k < mesh.nbt; k++) {
			Triangle tau = mesh.triangles[k];
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
		for (size_t k = 0; k < mesh.nbt; k++) {
			Triangle tau = mesh.triangles[k];
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
		for (size_t k = 0; k < mesh.nbt; k++) {
			Triangle tau = mesh.triangles[k];
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
//	for (size_t m = 0; m < mesh.nbv; m++) {
//		Vertex p = mesh.vertices[m];
//		double u = 0.0;
//		if (p.is_internal) {
//			size_t i = internal_ids[m];
//			u = x[i];
//		}
//		printf("%1.17le\t%1.17le\t%1.17le\n", p.pos.x, p.pos.y, u);
//	}
	for (size_t k = 0; k < mesh.nbt; k++) {
		Triangle tau = mesh.triangles[k];
		for (size_t v = 0; v < 4; v++) {
			double u = 0.0;
			Vertex p = *tau.vertices[v % 3];
			if (p.is_internal) {
				u = x[internal_ids[p.id]];
			}
			printf("%1.17le\t%1.17le\t%1.17le\n", p.pos.x, p.pos.y, u);
		}
		printf("\n\n");
	}

	free(x);
	free(b);
	free(a);
	free(a_data);
	free(internal_ids);
	free_mesh(mesh);
	return 0;
}

