// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

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
	char buffer[1024];
	while (fgets(buffer, sizeof(buffer), mesh_file) != NULL) {
		size_t id;
		double x, y; // Vertex (x, y)
		size_t n1, n2; // Edge (n1, n2)
		size_t e1, e2, e3; // Triangle (e1, e2, e3)
		if (sscanf(buffer, "Vertex %zu %lf %lf", &id, &x, &y) == 3) {
			printf("Vertex#%zu (%f, %f)\n", id, x, y);
		} else if (sscanf(buffer, "Edge %zu %zu %zu", &id, &n1, &n2) == 3) {
			printf("Edge#%zu (%zu, %zu)\n", id, n1, n2);
		} else if (sscanf(buffer, "Triangle %zu %zu %zu %zu", &id, &e1, &e2, &e3) == 4) {
			printf("Triangle#%zu (%zu, %zu, %zu)\n", id, e1, e2, e3);
		}
	}
	fclose(mesh_file);
	return 0;
}

// // // // // // // // // * キ ケ ン * // // // // // // // // // //
// // // // // malloc 関数の使用には十分注意すること！ // // // // //
// 事前に配列サイズが分からないときは，メモリ確保命令 malloc を使い
// 浮動小数点数型の配列を動的確保する．
// 使い終わったら必ず free() すること．
double * allocate_real_vector(size_t n) {
	double * vector = malloc(n * sizeof(double));
	if (vector == NULL) {
		fprintf(stderr, "Failed to allocate vector. ABORT.\n");
		exit(EXIT_FAILURE);
	}
	return vector;
}

