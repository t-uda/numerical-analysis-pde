// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀
//
// 正方形 [0, 1]^2 を合同な直角三角形で埋め尽くすように
// 一様に分割する．まず格子点 (i * dx, j * dy) を並べて
// 頂点とし，添え字同士の関係によって隣り合う頂点を辺とする．
// 最後に添え字間の関係で三頂点が繋がるところを三角形とする．
//
// 横の分割数 n, 縦の分割数 m としてメッシュ情報を出力する：
//
//     ./rectmesh.out n m > mesh.dat
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int main(int argc, char * argv[]) {
	size_t n = 10;
	if (argc >= 2) {
		n = strtoul(argv[1], NULL, 10);
	}
	size_t m = n;
	if (argc >= 3) {
		m = strtoul(argv[2], NULL, 10);
	}
	// 頂点数 (n+1)(m+1)
	printf("NbVertices %zu\n", (n + 1) * (m + 1));
	// 辺数 nm+n(m+1)+(n+1)m
	printf("NbEdges %zu\n", n * m + n * (m + 1) + (n + 1) * m);
	// 三角形数 2nm
	printf("NbTriangles %zu\n", 2 * n * m);
	printf("\n");
	double dx = 1.0 / n;
	double dy = 1.0 / m;
	// 各格子点に通し番号 v_id を振っていき，その座標と内部節店かどうかの真偽値を出力する．
	size_t v_id = 0;
	for (size_t i = 0; i <= n; ++i) {
		double x = i * dx;
		for (size_t j = 0; j <= m; ++j) {
			double y = j * dy;
			bool is_internal = (0 < i && i < n && 0 < j && j < m);
			printf("Vertex %zu %1.17lf %1.17lf %d\n", v_id, x, y, is_internal);
			v_id++; // v_id == i * (m + 1) + j
		}
	}
	printf("\n");
	// 各辺に通し番号 e_id を振っていき，二つの頂点番号を出力する．
	size_t e_id = 0;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			// 斜めの辺
			printf("Edge %zu %zu %zu\n", e_id, (i + 1) * (m + 1) + j, i * (m + 1) + j + 1);
			e_id++; // e_id == i * m + j
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j <= m; ++j) {
			// 横の辺
			printf("Edge %zu %zu %zu\n", e_id, i * (m + 1) + j, (i + 1) * (m + 1) + j);
			e_id++; // e_id == (n * m) + i * (m + 1) + j
		}
	}
	for (size_t i = 0; i <= n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			// 縦の辺
			printf("Edge %zu %zu %zu\n", e_id, i * (m + 1) + j, i * (m + 1) + j + 1);
			e_id++; // e_id == (n * m + n * (m + 1)) + i * m + j
		}
	}
	printf("\n");
	// 各三角形に通し番号 t_id を振っていき，三つの辺番号を出力する．
	size_t t_id = 0;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			size_t e0_id = i * m + j;
			size_t e1_id = n * m + i * (m + 1) + j;
			size_t e2_id = n * m + n * (m + 1) + i * m + j;
			// 左下の直角三角形
			printf("Triangle %zu %zu %zu %zu\n", t_id, e0_id, e1_id, e2_id);
			t_id++;
			// 向かい合う，直角を右上にもつ三角形
			printf("Triangle %zu %zu %zu %zu\n", t_id, e0_id, e1_id + 1, e2_id + m);
			t_id++;
		}
	}
	return 0;
}

