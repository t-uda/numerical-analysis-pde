// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀
//
// n×m の単純な三角形メッシュを生成する
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
	printf("NbVertices %zu\n", (n + 1) * (m + 1));
	printf("NbEdges %zu\n", n * m + n * (m + 1) + (n + 1) * m);
	printf("NbTriangles %zu\n", 2 * n * m);
	printf("\n");
	size_t v_id = 0;
	double dx = 1.0 / n;
	double dy = 1.0 / m;
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
	size_t e_id = 0;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			printf("Edge %zu %zu %zu\n", e_id, (i + 1) * (m + 1) + j, i * (m + 1) + j + 1);
			e_id++; // e_id == i * m + j
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j <= m; ++j) {
			printf("Edge %zu %zu %zu\n", e_id, i * (m + 1) + j, (i + 1) * (m + 1) + j);
			e_id++; // e_id == (n * m) + i * (m + 1) + j
		}
	}
	for (size_t i = 0; i <= n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			printf("Edge %zu %zu %zu\n", e_id, i * (m + 1) + j, i * (m + 1) + j + 1);
			e_id++; // e_id == (n * m + n * (m + 1)) + i * m + j
		}
	}
	printf("\n");
	size_t t_id = 0;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			size_t e0_id = i * m + j;
			size_t e1_id = n * m + i * (m + 1) + j;
			size_t e2_id = n * m + n * (m + 1) + i * m + j;
			printf("Triangle %zu %zu %zu %zu\n", t_id, e0_id, e1_id, e2_id);
			t_id++;
			printf("Triangle %zu %zu %zu %zu\n", t_id, e0_id, e1_id + 1, e2_id + m);
			t_id++;
		}
	}
	return 0;
}

