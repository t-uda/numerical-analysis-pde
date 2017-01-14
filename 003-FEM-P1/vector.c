
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vector.h"

// returns the subtraction a-b
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

