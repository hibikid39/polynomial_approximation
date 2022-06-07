#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define X_NUM 20
#define GROUP_NUM 5
#define START 0.0
#define END 1.0

#define N 11

/*
    [0, 1]の一様乱数からBox-Muller 法を用いてN(0, 0.05)正規乱数を生成
*/
double normrand() {
    double x1 = rand() / (double)RAND_MAX;  // [0,1]の一様分布
    double x2 = rand() / (double)RAND_MAX;  // [0,1]の一様分布
    double nr = sqrt(-2*log(x1))*cos(2*M_PI*x2);    // N(0, 1)の正規分布
    return nr*sqrt(0.05);    //　分散を0.05に調整
}

/*
    一方のベクトルからもう一方のベクトルを引く
    a, b : ベクトル
*/
void vec_diff(double *a, double *b) {
    for (int i = 0; i < N; i++) {
        b[i] -= a[i];
    }
}

/*
    ガウスの消去法によりN次連立方程式を解く
    m : 行列
    b : 定数項
*/
void Gauss(double **m, double *b) {
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            double coef = m[j][i] / m[i][i];
            double del[N];
            for (int k = 0; k < N; k++) {
                del[k] = m[i][k] * coef;
            }
            vec_diff(del, m[j]);
            b[j] -= b[i] * coef;
        }
    }

    for (int i = N -1; i >= 0; i--) {
        double x = 1. / m[i][i];
        m[i][i] *= x;
        b[i] *= x;
        for (int j = 0; j < i; j++) {
            b[j] -= b[i]*m[j][i];
            m[j][i] = 0;
        }
    }
}

/*
    多項式を計算する
    x : 引数
    w : 係数
*/
double calc_polynomial(double x, double *w) {
    double y = 0;
    for (int i = 0; i < N; i++) {
        y += w[N-1-i] * pow(x, (double)i);
    }
    return y;
}

/*
    行列を 出力
*/
void print_Matrix(double **m) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
    ベクトルを出力
*/
void print_Vector(double *v) {
    for (int i = 0; i < N; i++) {
        printf("%.2f ", v[i]);
    }
    printf("\n\n");
}

int main(int argc, char *argv[]) {
    int i, j, k, l;

    if (argc != 2) {
        printf("[Error] Usage g:generate data, a:approximate poly\n");
        return 0;
    }

    if (strcmp(argv[1], "g") == 0) // データ生成
    {
        FILE *fpw;
        fpw = fopen("data.txt", "w");
        if (fpw == NULL) {
            printf("[Error] could'nt open file");
            return 0;
        }

        double x[X_NUM], h[X_NUM], t[X_NUM];
        int seed = 1234567;
        srand(seed);

        for (j = 0; j < GROUP_NUM; j++) {
            for (i = 0; i < X_NUM; i++) {
                x[i] = (END - START) / (X_NUM - 1) * i; // データx
                h[i] = sin(2*M_PI*x[i]);                // データsin(2*PI*x)
                t[i] = h[i] + normrand();               // ノイズ追加
                fprintf(fpw, "%f %f\n", x[i], t[i]);
            }
            fprintf(fpw, "\n");
        }
    }
    else if (strcmp(argv[1], "a") == 0)    // 多項式近似を実行
    {
        // データ読み込み
        FILE *fpr;
        fpr = fopen("data.txt", "r");
        if (fpr == NULL) {
            printf("[Error] could'nt open file");
            return 0;
        }

        double x[X_NUM*GROUP_NUM], t[X_NUM*GROUP_NUM];
        for (j = 0; j < GROUP_NUM; j++) {
            for (i = 0; i < X_NUM; i++) {
                fscanf(fpr, "%lf%lf", &x[j*X_NUM + i], &t[j*X_NUM + i]);
            }
            fscanf(fpr, "\n");
        }

        // 交差確認法によりRMSの平均を求める
        double rms_sum_test = 0, rms_sum_train = 0;
        for (k = 0; k < GROUP_NUM; k++) {
            printf("--- Group %d ---\n", k+1);
            // 連立方程式を立てる （m, bを設定する）
            //double m[N][N];     // 行列M
            double **m = (double**)malloc(sizeof(double*)*N);
            for (i = 0; i < N; i++) m[i] = (double*)malloc(sizeof(double)*N);
            //double b[N];        // 定数項b
            double *b = (double*)malloc(sizeof(double)*N);;
            // mを設定
            double x_pow[2*N-1];
            for (i = 0; i < 2*N-1; i++) {
                double sum = 0;
                for(j = 0; j < X_NUM; j++) {
                    sum += pow(x[k*X_NUM + j], (double)i);
                }
                x_pow[i] = sum;
            }
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    m[i][j] = x_pow[2*(N-1)-(i+j)];
                }
            }
            printf("m = \n");
            print_Matrix(m);
            // bを設定
            for (i = 0; i < N; i++) {
                double sum = 0;
                for (j = 0; j< X_NUM; j++) {
                    sum += pow(x[k*X_NUM + j], (double)(N-1-i)) * t[k*X_NUM + j];
                }
                b[i] = sum;
            }
            printf("b = \n");
            print_Vector(b);

            // 連立方程式をガウスの消去法により解く
            Gauss(m, b);
            printf("w = \n");
            print_Vector(b);

            // テストデータの二乗誤差
            for (l = 0; l < GROUP_NUM; l++) {
                if (l != k) {
                    double y;
                    for (i = 0; i < X_NUM; i++) {
                        y = calc_polynomial(x[l*X_NUM + i], b);
                        rms_sum_test += (t[l*X_NUM + i] - y)*(t[l*X_NUM + i] - y);
                        //printf("%f, %f, %f\n", x[l*X_NUM + i], t[l*X_NUM + i], calc_polynomial(x[l*X_NUM + i], &b));
                    }
                }
            }

            // 学習データの平均二乗誤差
            double y;
            for (i = 0; i < X_NUM; i++) {
                y = calc_polynomial(x[k*X_NUM + i], b);
                rms_sum_train += (t[k*X_NUM + i] - y)*(t[k*X_NUM + i] - y);
                //printf("%f, %f, %f\n", x[k*X_NUM + i], t[k*X_NUM + i], calc_polynomial(x[l*X_NUM + i], &b));
            }
        }


        printf("rms_train, rms_test = %f, %f\n", sqrt(rms_sum_train/X_NUM), sqrt(rms_sum_test/((GROUP_NUM-1)*X_NUM)));

    }

    printf("End.\n");
    return 0;
}
