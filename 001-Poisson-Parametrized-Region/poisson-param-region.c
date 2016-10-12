// 例1.2 境界が媒介変数表示された領域でのポアソン方程式の差分解法．
//
// ここでは単位円盤領域において厳密解が
//
//     u = x^2 - y^2
//
// となるラプラス方程式を解いている．つまり，ポアソン方程式の
// 右辺はここでは 0 とし，境界条件として x^2 - y^2 を与えている．
//
// プログラムは右辺の関数 f およびディリクレデータが任意の場合に
// 対応しているので，自由に変えて試してみよ．また，凸領域であって
// 上下左右の局所座標表示が陽に与えられる場合にも対応しているので，
// 例えば楕円型領域などに変えて試してみよ．
//
// プログラムは make を使うか，gcc を使ってコンパイルせよ．
// このプログラムは標準出力に計算結果を出力する．
// 可視化するためには，一度ファイルに出力するなどして，
// 例えば gnuplot で splot 'output.dat' を試せ．
//
//     $ make
//     $ ./poisson-param-region.out > output.dat
//     $ gnuplot -e "splot 'output.dat'"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// 矩形領域 [X0, X0 + L] \times [Y0, Y0 + L] に格子点を作る
// 格子の基準点
#define X0 (-1.0)
#define Y0 (-1.0)
// 格子の縦および横の長さ
#define L (2.0)
// L あたりの分割数
#define N (128)
// 分割幅
#define h (L / N)

// ヤコビ法の反復終了の閾値
#define TOLERANCE (1.0e-13)
// ヤコビ法の最大反復回数
#define MAXITER (20000)
// 出力時に間引くデータ間隔
#define CULL (4)

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

// 境界のパラメータ表示
// 凸領域でない場合，上下左右それぞれの表示が一意に定まるとは限らないことに注意
double y_top(double x) {
	return sqrt(1.0 - x * x);
}
double y_bottom(double x) {
	return -sqrt(1.0 - x * x);
}
double x_right(double y) {
	return sqrt(1.0 - y * y);
}
double x_left(double y) {
	return -sqrt(1.0 - y * y);
}

// 内部節点かどうか判定する．
// 境界は false を返す．
bool is_internal(int i, int j) {
	double x = X0 + i * h;
	double y = Y0 + j * h;
	if (y >= y_top(x)) return false;
	if (y <= y_bottom(x)) return false;
	if (x >= x_right(y)) return false;
	if (x <= x_left(y)) return false;
	return true;
}

// 内部節点 P_ij における lambda_X (X=B,C,D,E) の値を返す．
// 凸領域で上下左右の媒介変数表示が与えられているため以下のように
// 簡単に求まるが，一般にこの方法で求められるとは限らない．
// （つまり，より丁寧な場合分けが必要である．）
// 右側
double lambda_B(int i, int j) {
	if (is_internal(i + 1, j)) return 1.0;
	return x_right(Y0 + j * h) - (X0 + i * h);
}
// 左側
double lambda_C(int i, int j) {
	if (is_internal(i - 1, j)) return 1.0;
	return (X0 + i * h) - x_left(Y0 + j * h);
}
// 上側
double lambda_D(int i, int j) {
	if (is_internal(i, j + 1)) return 1.0;
	return y_top(X0 + i * h) - (Y0 + j * h);
}
// 下側
double lambda_E(int i, int j) {
	if (is_internal(i, j - 1)) return 1.0;
	return (Y0 + j * h) - y_bottom(X0 + i * h);
}

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
	// // // // // malloc 関数の使用には十分注意すること！ // // // // //
	// 事前に未知数の数 n が分からないため，メモリ確保命令 malloc を使い
	// 浮動小数点数型の配列を動的確保する．
	double * uh = malloc(n * sizeof(double));
	double * uh_new = malloc(n * sizeof(double));
	double * Fh = malloc(n * sizeof(double));
	if (uh == NULL || uh_new == NULL || Fh == NULL) {
		fprintf(stderr, "Failed to allocate vectors. ABORT.\n");
		free(uh);
		free(uh_new);
		free(Fh);
		exit(EXIT_FAILURE);
	}

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
					double y = Y0 + j * h;
					Fh[m] += Dirichlet_data(x_right(y), y) / (lambda_B(i, j) * h * h);
				}
				if (!is_internal(i - 1, j)) { // 左側
					double y = Y0 + j * h;
					Fh[m] += Dirichlet_data(x_left(y), y) / (lambda_C(i, j) * h * h);
				}
				if (!is_internal(i, j + 1)) { // 上側
					double x = X0 + i * h;
					Fh[m] += Dirichlet_data(x, y_top(x)) / (lambda_D(i, j) * h * h);
				}
				if (!is_internal(i, j - 1)) { // 下j側
					double x = X0 + i * h;
					Fh[m] += Dirichlet_data(x, y_bottom(x)) / (lambda_E(i, j) * h * h);
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
				// debug
				if (fabs(X0 + i * h) < h / 8 && fabs(Y0 + j * h) < h / 8) {
					fprintf(stderr, "u(0, 0) = %e\n", uh[m]);
				}
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

