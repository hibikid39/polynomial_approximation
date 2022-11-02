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
    [0, 1]�̈�l��������Box-Muller �@��p����N(0, 0.05)���K�����𐶐�
*/
double normrand() {
    double x1 = rand() / (double)RAND_MAX;  // [0,1]�̈�l���z
    double x2 = rand() / (double)RAND_MAX;  // [0,1]�̈�l���z
    double nr = sqrt(-2*log(x1))*cos(2*M_PI*x2);    // N(0, 1)�̐��K���z
    return nr*sqrt(0.05);    //�@���U��0.05�ɒ���
}

/*
    ����̃x�N�g�������������̃x�N�g��������
    a, b : �x�N�g��
*/
void vec_diff(double *a, double *b) {
    for (int i = 0; i < N; i++) {
        b[i] -= a[i];
    }
}

/*
    �K�E�X�̏����@�ɂ��N���A��������������
    m : �s��
    b : �萔��
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
    ���������v�Z����
    x : ����
    w : �W��
*/
double calc_polynomial(double x, double *w) {
    double y = 0;
    for (int i = 0; i < N; i++) {
        y += w[N-1-i] * pow(x, (double)i);
    }
    return y;
}

/*
    �s��� �o��
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
    �x�N�g�����o��
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

    if (strcmp(argv[1], "g") == 0) // �f�[�^����
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
                x[i] = (END - START) / (X_NUM - 1) * i; // �f�[�^x
                h[i] = sin(2*M_PI*x[i]);                // �f�[�^sin(2*PI*x)
                t[i] = h[i] + normrand();               // �m�C�Y�ǉ�
                fprintf(fpw, "%f %f\n", x[i], t[i]);
            }
            fprintf(fpw, "\n");
        }
    }
    else if (strcmp(argv[1], "a") == 0)    // �������ߎ������s
    {
        // �f�[�^�ǂݍ���
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

        // �����m�F�@�ɂ��RMS�̕��ς����߂�
        double rms_sum_test = 0, rms_sum_train = 0;
        for (k = 0; k < GROUP_NUM; k++) {
            printf("--- Group %d ---\n", k+1);
            // �A���������𗧂Ă� �im, b��ݒ肷��j
            //double m[N][N];     // �s��M
            double **m = (double**)malloc(sizeof(double*)*N);
            for (i = 0; i < N; i++) m[i] = (double*)malloc(sizeof(double)*N);
            //double b[N];        // �萔��b
            double *b = (double*)malloc(sizeof(double)*N);;
            // m��ݒ�
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
            // b��ݒ�
            for (i = 0; i < N; i++) {
                double sum = 0;
                for (j = 0; j< X_NUM; j++) {
                    sum += pow(x[k*X_NUM + j], (double)(N-1-i)) * t[k*X_NUM + j];
                }
                b[i] = sum;
            }
            printf("b = \n");
            print_Vector(b);

            // �A�����������K�E�X�̏����@�ɂ�����
            Gauss(m, b);
            printf("w = \n");
            print_Vector(b);

            // �e�X�g�f�[�^�̓��덷
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

            // �w�K�f�[�^�̕��ϓ��덷
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
