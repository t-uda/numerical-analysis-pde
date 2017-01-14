// * 京都大学 2016 年度後期「数値解析」
// 担当教員：Karel Svadlenka
// プログラム作成者：宇田 智紀
//
// n×m の単純な三角形メッシュを生成する
//     ./genmesh n m > mesh.dat
//

#include <stdlib.h>

#include "mesh.h"

double g(double x, double y) {
	return x * x + y * y - 1.0;
}

int main(int argc, char * argv[]) {
	return 0;
}

