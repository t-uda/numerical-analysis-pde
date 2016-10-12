// 例1.1 正方形領域におけるポアソン方程式の差分解法．
//
// ここではディリクレ境界 0 で，厳密解が
//
//     u = \sin(\pi x) \sin(\pi y).
//
// となるポアソン方程式を解いている．
// プログラムは make を使うか，gcc を使ってコンパイルせよ．
// このプログラムは標準出力に計算結果を出力する．
// 可視化するためには，一度ファイルに出力するなどして，
// 例えば gnuplot で splot 'output.dat' を試せ．
//
//     $ make
//     $ ./poisson-rect.out > output.dat
//     $ gnuplot -e "splot 'output.dat'"

#include <stdio.h>
#include <math.h>

// 領域の縦横分割数
#define N (200)
// 領域の縦横分割幅
#define h (1.0 / N)
// 未知変数の数（内部節点の数）
#define n ((N-1)*(N-1))

// ヤコビ法の反復終了の閾値
#define TOLERANCE (1.0e-10)
// ヤコビ法の最大反復回数
#define MAXITER (25000)
// 出力時に間引くデータ間隔
#define CULL (5)

double uh[n];
double uh_new[n];
double Fh[n];

// 内部節点 (ih, jh) (0 < i, j < N) に対して 0 <= m < n を返す
int get_index(int i, int j) {
	return (N - 1) * (i - 1) + (j - 1);
}

// ポアソン方程式の右辺
double f(double x, double y) {
	return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

int main(int argc, char * argv[]) {
	// 右辺ベクトルの初期化
	for (int i = 1; i < N; ++i) {
		for (int j = 1; j < N; ++j) {
			// ディリクレ zero の場合の右辺ベクトルは単に f である．
			// 講義でも述べられていたように，ディリクレ nonzero の場合は
			// 境界条件に由来する項を右辺ベクトルに足し上げよ．
			int m = get_index(i, j);
			Fh[m] = f(i * h, j * h);
		}
	}

	// ヤコビ法による反復計算で方程式 Ah * uh = Fh を解く
	// 反復の初期値は D^{-1} * Fh で与える
	for (int m = 0; m < n; ++m) {
		uh[m] = Fh[m] / (4 * N * N);
	}
	for (int k = 0; k < MAXITER; ++k) {
		// 疎行列 H = -D^{-1} * R は非零要素 1/4 を非対角に高々 4 つもつ
		for (int i = 1; i < N; ++i) {
			for (int j = 1; j < N; ++j) {
				int m = get_index(i, j);
				uh_new[m] = Fh[m] / (4 * N * N); // D^{-1} * Fh
				if (i > 1) {
					uh_new[m] += uh[get_index(i - 1, j)] / 4.0;
				}
				if (i < N - 1) {
					uh_new[m] += uh[get_index(i + 1, j)] / 4.0;
				}
				if (j > 1) {
					uh_new[m] += uh[get_index(i, j - 1)] / 4.0;
				}
				if (j < N - 1) {
					uh_new[m] += uh[get_index(i, j + 1)] / 4.0;
				}
			}
		}

		// 未知ベクトルを更新する
		double uh_norm = 0.0;
		double res = 0.0;
		for (int i = 1; i < N; ++i) {
			for (int j = 1; j < N; ++j) {
				int m = get_index(i, j);
				res += (uh[m] - uh_new[m]) * (uh[m] - uh_new[m]);
				uh_norm += uh[m] * uh[m];
				uh[m] = uh_new[m];
			}
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
			if (i == 0 || i == N || j == 0 || j == N) {
				printf("%f %f %f\n", i * h, j * h, 0.0);
			} else {
				int m = get_index(i, j);
				printf("%f %f %f\n", i * h, j * h, uh[m]);
			}
		}
	}
	return 0;
}

