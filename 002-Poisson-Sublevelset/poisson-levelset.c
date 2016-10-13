// 例1.2' 劣位集合で書かれる領域上のポアソン方程式の差分解法．
//
// このプログラムは，領域が劣位集合，つまり {g(x, y) < 0} の
// 形で書かれる一般の場合に対応している．それ以外の点は
// 前のプログラム 001-Poisson-Parametrized-Region と同様である．
// 特に，非凸領域にも対応している．
//
// プログラムは make を使うか，gcc を使ってコンパイルせよ．
// このプログラムは標準出力に計算結果を出力する．
// 可視化するためには，一度ファイルに出力するなどして，
// 例えば gnuplot で splot 'output.dat' を試せ．
//
//     $ make
//     $ ./poisson-levelset.out > output.dat
//     $ gnuplot -e "splot 'output.dat'"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// 矩形領域 [X0, X0 + L] \times [Y0, Y0 + L] に格子点を作る
// 格子の基準点
#define X0 (-1.5)
#define Y0 (-1.5)
// 格子の縦および横の長さ
#define L (3.0)
// L あたりの分割数
#define N (256)
// 分割幅
#define h (L / N)

// ヤコビ法の反復終了の閾値
#define TOLERANCE (1.0e-12)
// ヤコビ法の最大反復回数
#define MAXITER (10000)
// 出力時に間引くデータ間隔
#define CULL (2)

// 内部節点に連番をふるための配列
int indices[N + 1][N + 1];

// ポアソン方程式の右辺
double f(double x, double y) {
	return 0.0;
}

// 境界条件
double Dirichlet_data(double x, double y) {
	return x * x - y * y;
}

// 境界を定める陰関数；領域内部は劣位集合 {g(x, y) < 0} で与えられるものとする
double g(double x, double y) {
	return (x * x + y * y) * (2.0 + cos(M_PI * x) + cos(M_PI * y)) - 1.15;
}
// g の偏微分
double dg_dx(double x, double y) {
	return 2.0 * x * (2.0 + cos(M_PI * x) + cos(M_PI * y)) - M_PI * sin(M_PI * x) * (x * x + y * y);
}
double dg_dy(double x, double y) {
	return 2.0 * y * (2.0 + cos(M_PI * x) + cos(M_PI * y)) - M_PI * sin(M_PI * y) * (x * x + y * y);
}

// g(x, y0) = 0 なる x を Newton 法で求める
double compute_boundary_x(double x, double y0) {
	for (size_t k = 0; k < 10; ++k) {
		double dx = -g(x, y0) / dg_dx(x, y0);
		if (fabs(dx) < 1.0e-15) break;
		x += dx;
	}
	return x;
}

// g(x0, y) = 0 なる y を Newton 法で求める
double compute_boundary_y(double x0, double y) {
	for (size_t k = 0; k < 10; ++k) {
		double dy = -g(x0, y) / dg_dy(x0, y);
		if (fabs(dy) < 1.0e-15) break;
		y += dy;
	}
	return y;
}

// 内部節点かどうか判定する．
// 境界は false を返す．
bool is_internal(int i, int j) {
	double x = X0 + i * h;
	double y = Y0 + j * h;
	return (g(x, y) < 0);
}

// 内部節点 P_ij における lambda_X (X=B,C,D,E) の値を返す．
// 右側
double lambda_B(int i, int j) {
	if (is_internal(i + 1, j)) return 1.0;
	return compute_boundary_x(X0 + i * h, Y0 + j * h) - (X0 + i * h);
}
// 左側
double lambda_C(int i, int j) {
	if (is_internal(i - 1, j)) return 1.0;
	return (X0 + i * h) - compute_boundary_x(X0 + i * h, Y0 + j * h);
}
// 上側
double lambda_D(int i, int j) {
	if (is_internal(i, j + 1)) return 1.0;
	return compute_boundary_y(X0 + i * h, Y0 + j * h) - (Y0 + j * h);
}
// 下側
double lambda_E(int i, int j) {
	if (is_internal(i, j - 1)) return 1.0;
	return (Y0 + j * h) - compute_boundary_y(X0 + i * h, Y0 + j * h);
}

// 浮動小数点数型の配列を動的確保する関数．後で定義する．
double * allocate_real_vector(size_t);

int main(int argc, char * argv[]) {
	// 内部節点の数 n を数える（未知数の数）
	int n = 0;
	for (int i = 0; i <= N; ++i) {
		for (int j = 0; j <= N; ++j) {
			if (is_internal(i, j)) {
				// 節点に順番に番号を振って格納しておく
				indices[i][j] = n;
				++n;
			} else {
				indices[i][j] = -1;
			}
		}
	}

	// // // // // // // // // * キ ケ ン * // // // // // // // // // //
	// 事前に未知数の数 n が分からないので配列は動的確保する．
	// 動的確保した配列は free() を使って解放するのを忘れないこと．
	double * uh = allocate_real_vector(n);
	double * uh_new = allocate_real_vector(n);
	double * Fh = allocate_real_vector(n);

	// 右辺ベクトルを初期化する
	for (int i = 0; i <= N; ++i) {
		for (int j = 0; j <= N; ++j) {
			if (is_internal(i, j)) {
				int m = indices[i][j];
				double x = X0 + i * h;
				double y = Y0 + j * h;
				// 正則な節点の場合，C = 1 となることに注意
				double C = (lambda_B(i, j) + lambda_C(i, j)) * (lambda_D(i, j) + lambda_E(i, j)) / 4.0;
				Fh[m] = C * f(x, y);
				// 境界条件から既知の数を右辺ベクトルに足し上げる．
				if (!is_internal(i + 1, j)) { // 右側
					Fh[m] += Dirichlet_data(compute_boundary_x(x, y), y) / (lambda_B(i, j) * h * h);
				}
				if (!is_internal(i - 1, j)) { // 左側
					Fh[m] += Dirichlet_data(compute_boundary_x(x, y), y) / (lambda_C(i, j) * h * h);
				}
				if (!is_internal(i, j + 1)) { // 上側
					Fh[m] += Dirichlet_data(x, compute_boundary_y(x, y)) / (lambda_D(i, j) * h * h);
				}
				if (!is_internal(i, j - 1)) { // 下側
					Fh[m] += Dirichlet_data(x, compute_boundary_y(x, y)) / (lambda_E(i, j) * h * h);
				}
			}
		}
	}

	// ヤコビ法による反復計算で方程式 Ah * uh = Fh を解く
	// 反復の初期値は D^{-1} * Fh に近いものでよい
	for (int m = 0; m < n; ++m) {
		uh[m] = Fh[m] * h * h / 4.0;
	}

	for (int k = 0; k < MAXITER; ++k) {
		// 疎行列 H = -D^{-1} * R は非零要素を非対角に高々 4 つもつ
		for (int i = 0; i <= N; ++i) {
			for (int j = 0; j <= N; ++j) {
				if (is_internal(i, j)) {
					int m = indices[i][j];
					double L_B = 1.0 / lambda_B(i, j);
					double L_C = 1.0 / lambda_C(i, j);
					double L_D = 1.0 / lambda_D(i, j);
					double L_E = 1.0 / lambda_E(i, j);
					double a_ii = L_B + L_C + L_D + L_E;
					uh_new[m] = Fh[m] * h * h / a_ii; // D^{-1} * Fh
					if (is_internal(i + 1, j)) { // 右側
						uh_new[m] += uh[indices[i + 1][j]] * L_B / a_ii;
					}
					if (is_internal(i - 1, j)) { // 左側
						uh_new[m] += uh[indices[i - 1][j]] * L_C / a_ii;
					}
					if (is_internal(i, j + 1)) { // 上側
						uh_new[m] += uh[indices[i][j + 1]] * L_D / a_ii;
					}
					if (is_internal(i, j - 1)) { // 下側
						uh_new[m] += uh[indices[i][j - 1]] * L_E / a_ii;
					}
				}
			}
		}

		// 未知ベクトルを更新する
		double uh_norm = 0.0;
		double res = 0.0;
		for (int m = 0; m < n; ++m) {
			res += (uh[m] - uh_new[m]) * (uh[m] - uh_new[m]);
			uh_norm += uh[m] * uh[m];
			uh[m] = uh_new[m];
		}
		res /= uh_norm;
		// 更新分の二乗和を見て反復計算を続けるか判断する
		if (res < TOLERANCE || k == MAXITER - 1) {
			// 標準エラー出力にログを書き出す．
			// このようにデバッグ用のログは標準エラー出力に出しておけば，
			// 標準出力にでた数値結果の方を汚さずにファイルに保存できる．
			fprintf(stderr, "[%d] res = %e\n", k, res);
			break;
		}
	}

	// 出力する
	for (int i = 0; i <= N; i += CULL) {
		for (int j = 0; j <= N; j += CULL) {
			if (is_internal(i, j)) {
				int m = indices[i][j];
				printf("%f %f %f\n", X0 + i * h, Y0 + j * h, uh[m]);
			}
		}
	}

	// // // // // // // // // // * 重 要 * // // // // // // // // // //
	// // // // // malloc で確保したメモリは必ず解放する！ // // // // //
	free(uh);
	free(uh_new);
	free(Fh);
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

