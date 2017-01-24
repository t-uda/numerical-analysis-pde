// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀
//
// 劣位集合 { (x, y) | g(x, y) < 0 } の形で与えられた領域を，
// 正三角形で埋め尽くすように分割する．プログラムではまず，
// 劣位集合内に正三角形を一定の規則に従って並べメッシュを作る．
// 次に，境界に近い位置にある頂点を境界上に移動（近似）する．
//
// 横の分割数 n, 縦の分割数 m としてメッシュ情報を出力する：
//
//     ./genmesh.out n m > mesh.dat
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"
#include "mesh.h"

// 頂点はだいたい矩形 (X0, X0 + L) x (Y0, Y0 + L') 内に配置する．
#define L (2.10)
#define X0 (-1.05)
#define Y0 (-1.05)

// 領域を劣位集合 { (x, y) | g(x, y) < 0 } によって与える
double g(double x, double y) {
	return x * x + y * y - 1.0;
}

// g の勾配
Vector2D grad_g(double x, double y) {
	Vector2D grad = {2.0 * x, 2.0 * y};
	return grad;
}

// p から dir 方向にもっとも近い境界点 q を Newton 法で求める；
// Find a real value s such that g(p + s * dir) = 0.
Vector2D compute_boundary(Vector2D p, Vector2D dir) {
	Vector2D q;
	double s = 0.0;
	for (size_t iter = 0; iter < 10; iter++) {
		// q := p + s * dir
		q = vec2D_add(p, vec2D_mul(s, dir));
		// s := s - g(q) / <dir, grad g(q)>
		double d = g(q.x, q.y) / vec2D_prod(dir, grad_g(q.x, q.y));
		if (fabs(d) < 1.0e-15) break;
		s -= d;
	}
	return q;
}

int main(int argc, char * argv[]) {
	size_t n = (argc >= 2) ?  strtoul(argv[1], NULL, 10) : 32;
	size_t m = (argc >= 3) ?  strtoul(argv[2], NULL, 10) : ceil(n * 2.0 / sqrt(3.0));
	double h = L / n;

	Mesh mesh;

	// 頂点は高々 (n-1)m 個
	mesh.vertices = malloc((n - 1) * m * sizeof(Vertex));
	check_allocated_or_abort(mesh.vertices);

	// 頂点へのポインタを格納する二次元配列（添え字番号は下記参照）
	// 全ての要素は NULL で初期化されている
	Vertex ** table_data = calloc((n + 1) * (m + 1), sizeof(Vertex *));
	check_allocated_or_abort(table_data);
	Vertex *** table = malloc((n + 1) * sizeof(Vertex **));
	check_allocated_or_abort(table);
	for (size_t k = 0; k <= n; k++) {
		table[k] = &table_data[k * (m + 1)];
	}

	{
		// (1, 0), (-1/2, sqrt(3)/2) を基底とする斜交座標 (h * i, h * j) を用いる．
		// 頂点はあとからこの座標を用いて table[k][j] で参照できるようにする．
		size_t id = 0;
		for (size_t j = 0; j < m; j++) {
			for (size_t k = 1; k < n; k++) {
				size_t i = j / 2 + k;
				double x = X0 + h * (i - 0.5 * j);
				double y = Y0 + h * (sqrt(0.75) * j);
				if (g(x, y) < 0.0) {
					mesh.vertices[id].id = id;
					mesh.vertices[id].pos.x = x;
					mesh.vertices[id].pos.y = y;
					mesh.vertices[id].is_internal = (g(x, y) < 0.0);
					table[k][j] = &mesh.vertices[id];
					id++;
				}
			}
		}
		mesh.nbv = id;
	}

	// 辺は高々 3(n-1)m 本
	mesh.edges = malloc(3 * (n - 1) * m * sizeof(Edge));
	check_allocated_or_abort(mesh.edges);

	{
		size_t id = 0;
		// 上と同様にループを回し，テーブルを参照して隣り合わせに頂点があれば辺を登録する
		for (size_t j = 0; j < m; j++) {
			for (size_t k = 1; k < n; k++) {
				if (table[k][j] == NULL) continue;
				if (table[k+1][j] != NULL) { // 右隣
					mesh.edges[id].p = table[k][j];
					mesh.edges[id].q = table[k+1][j];
					id++;
				}
				if (table[k-(j%2)][j+1] != NULL) { // 左上（jの偶奇で添え字がずれる）
					mesh.edges[id].p = table[k][j];
					mesh.edges[id].q = table[k-(j%2)][j+1];
					id++;
				}
				if (table[k+1-(j%2)][j+1] != NULL) { // 右上（jの偶奇で添え字がずれる）
					mesh.edges[id].p = table[k][j];
					mesh.edges[id].q = table[k+1-(j%2)][j+1];
					id++;
				}
			}
		}
		mesh.nbe = id;
	}

	// 三角形は高々 2(n-1)m 個
	mesh.triangles = malloc(2 * (n - 1) * m * sizeof(Triangle));
	check_allocated_or_abort(mesh.triangles);

	{
		size_t id = 0;
		// 上と同様にループを回し，テーブルを参照して三角形を形成していれば登録する
		for (size_t j = 0; j < m; j++) {
			for (size_t k = 1; k < n; k++) {
				Vertex * p = table[k][j];
				if (p == NULL) continue;
				Vertex * q = table[k+1][j]; // 右隣
				Vertex * r = table[k-(j%2)][j+1]; // 左上
				Vertex * s = table[k+1-(j%2)][j+1]; // 右上
				if (q != NULL && s != NULL) { // 三角形pqs
					Triangle * tau = &mesh.triangles[id];
					tau->vertices[0] = p;
					tau->vertices[1] = q;
					tau->vertices[2] = s;
					tau->edges[0] = find_edge(mesh.nbe, mesh.edges, p, q);
					tau->edges[1] = find_edge(mesh.nbe, mesh.edges, q, s);
					tau->edges[2] = find_edge(mesh.nbe, mesh.edges, s, p);
					id++;
				}
				if (r != NULL && s != NULL) { // 三角形prs
					Triangle * tau = &mesh.triangles[id];
					tau->vertices[0] = p;
					tau->vertices[1] = r;
					tau->vertices[2] = s;
					tau->edges[0] = find_edge(mesh.nbe, mesh.edges, p, r);
					tau->edges[1] = find_edge(mesh.nbe, mesh.edges, r, s);
					tau->edges[2] = find_edge(mesh.nbe, mesh.edges, s, p);
					id++;
				}
			}
		}
		mesh.nbt = id;
	}

	{
		// 今のままでは劣位集合の内部に多角形分割を作っただけなので，
		// 多角形の境界辺を正しい境界に近づくよう移動する．
		for (size_t e_id = 0, nbe = mesh.nbe; e_id < nbe; e_id++) {
			Edge * e = &mesh.edges[e_id];
			Triangle * adj_tau = NULL;
			size_t count = 0;
			// 辺 e をもつ三角形の個数を数え上げる
			for (size_t t_id = 0; t_id < mesh.nbt; t_id++) {
				Triangle * tau = &mesh.triangles[t_id];
				if (e == tau->edges[0] || e == tau->edges[1] || e == tau->edges[2]) {
					count++;
					adj_tau = tau;
				}
			}
			// 1 個の三角形にしか面していない辺は境界にある
			if (count == 1) {
				Vector2D p = e->p->pos, q = e->q->pos;
				e->p->pos = compute_boundary(p, grad_g(p.x, p.y));
				e->p->is_internal = false;
				e->q->pos = compute_boundary(q, grad_g(q.x, q.y));
				e->q->is_internal = false;
			}
		}
	}

	// mesh を標準出力に書き出す
	write_mesh(mesh, stdout);

	free(table_data);
	free(table);
	free_mesh(mesh);
	return 0;
}

