#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <fstream>

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define Nnot	10
#define NumMacro 1000
#define NumSys	1
#define NumConstraint	1
#define NumThreshold	1

// inputs for Generate R(0,1) by L'ecuyer (1997)
#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0


double MRG32k3a(void);  //Generate R(0,1) by L'ecuyer (1997)
// choices of seeds for Generate R(0,1) by L'ecuyer (1997)
double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;
//double  s10 = 43, s11 =54, s12 =65, s20 =43, s21 =54, s22 =65;
//double  s10 = 4321111, s11 =1115432, s12 =1116543, s20 =4321111, s21 =1115432, s22 =6543111;
//double  s10 = 43221, s11 =54332, s12 =65443, s20 =43321, s21 =54532, s22 =61543;
//double  s10 = 1010, s11 =10, s12 =101, s20 =2001, s21 = 202, s22 = 202;

double minfn(double x1, double x2);
double maxfn(double x, double y);

//double normal(double rmean, double rvar);
double configuration(void);
//double generate_multiNormal(int numConstraint, int case_index);
double generate_Bernoulli(int numConstraint);
//int read_chol_matrix(void);

double mean_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];
double q[NumThreshold][NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];

int main()
{
    double total_obs;   //총 obs 갯수
    double final_cd;    //correct decision 인지 아닌지
    double correct_decision;    // 최종 correct decision 갯수
    double overall_obs;     // 최종 obs 갯수

    correct_decision = 0;
    overall_obs = 0;

    double alpha = 0.05;
    double k = NumSys;
    double s = NumConstraint;
    double n0 = Nnot;
    double beta[NumSys];
    double eta[NumConstraint];

    beta[0] = (1 - pow(1 - alpha, 1 / k)) / s;
    printf("beta: %.10f\n", beta[0]);
    for (int j = 0; j < NumConstraint; j++) {
        eta[j] = (pow(2 * beta[0], -1 * (2 / (n0 - 1))) - 1) / 2;
        printf("eta: %.10f\n", eta[j]);
    }

    // >=2 threshold
    /*for (int j = 0; j < NumConstraint; j++) {
        eta[j] = (pow(2 * beta[0], -1 * (2 / (n0 - 1))) - 1) / 2;
        printf("eta: %.10f\n", eta[j]);
    }*/

    for (int l = 0; l < NumMacro; l++) {    // iteration 시작. 

        configuration();    //최초설정
        total_obs = 0;      //우선 obs 갯수 0으로 시작
        final_cd = 1;       //일단 cd를 1로 시작

        double num_obs[NumSys][NumConstraint];
        double R[NumSys][NumConstraint];
        double Sil2[NumSys][NumConstraint];

        for (int i = 0; i < NumSys; i++) {
            for (int j = 0; j < NumConstraint; j++) {
                Sil2[i][j] = 0;     //l번째 iter에서 i 시스템의 표본분산
                num_obs[i][j] = 0;      //i번째 시스템의 j 번째 제약에서 갯수
            }
        }

        for (int i = 0; i < NumSys; i++) {

            // generate initial samples
            double sumY[NumConstraint];     //obs 값들의 합
            double sum_squareY[NumConstraint];      //obs값들의 sum square
            //double sumYq[NumConstraint];    //obs 값 - q의 값들의 합

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] = 0;
                sum_squareY[j] = 0;
            }

            for (int n = 0; n < Nnot; n++) {    //Nnot 만큼 최초 obs 생성
                generate_Bernoulli(NumConstraint);  // obs 하나 만듦 (베르누이 분포에서 평균)
                total_obs += 1;     // constraint와 상관없이 시스템의 총 obs 갯수 모든 시스템에 대해 같이 증가함.

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    sum_squareY[j] += observations[i][j] * observations[i][j];
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    //mean_value[i][j] = (sumY[j] / Nnot);
                }


                for (int j = 0; j < NumConstraint; j++) surviveThreshold[j] = NumThreshold;     //threshold 개수로 survive한거 정의 하고 시작

            }

            // find continuation region. 일단 Nnot까지 한거에 대해서 구함
            for (int j = 0; j < NumConstraint; j++) {
                Sil2[i][j] = (sum_squareY[j] / (Nnot - 1)) - (sumY[j] / (Nnot - 1)) * (sumY[j] / Nnot);
                R[i][j] = maxfn(0, (Nnot - 1) * Sil2[i][j] * (eta[j]) / epsilon[j] - epsilon[j] * num_obs[i][j] / 2);
                //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];    //평행선 continuation
            }

            while (surviveConstraint != 0) {

                for (int j = 0; j < NumConstraint; j++) {   //모든 const에 대해서 아래 검사

                    if (ON[i][j] == 1) {    //만약 그 시스템에서 그 constraint가 아직 검사가 안됐다면

                        for (int d = 0; d < NumThreshold; d++) {    //threshold까지 내려가서 검사

                            if (ON_l[i][j][d] == 1) {   //검사 안됐었으면

                                if ((sumY[j] + R[i][j]) / num_obs[i][j] <= q[d][j]) {     // feasible 조건
                                    Z[i][j][d] = 1;     // i 시스템의 j 제약의 d번째 threshold에 대해 feasible
                                    ON_l[i][j][d] = 0;  // 해당 threshold를 검사한거로 변경
                                    surviveThreshold[j] -= 1;
                                }

                                if ((sumY[j] - R[i][j]) / num_obs[i][j] >= q[d][j]) {   // infeasible 조건
                                    Z[i][j][d] = 0;     // i 시스템의 j 제약의 d번째 threshold에 대해 infeasible
                                    ON_l[i][j][d] = 0;      // 해당 threshold를 검사한거로 변경
                                    surviveThreshold[j] -= 1;
                                }
                            }   //아무거에도 해당 안됐었으면 그냥 다음으로 넘어가고 obs 하나 더 얻는다 
                        }

                        if (surviveThreshold[j] == 0) {     //더이상 검사안한 threshold가 없으면 해당 constraint도 검사한거로 변경
                            ON[i][j] = 0;
                            surviveConstraint -= 1;
                        }

                    }

                }

                if (surviveConstraint == 0) break;  //검사 안한 제약이 없으면 종료

                generate_Bernoulli(NumConstraint);      //만약 위 조건을 만족 안했었다면 obs 하나 더 생성
                total_obs += 1;     //시스템의 total obs도 하나 늘리기

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    num_obs[i][j] += 1;
                    R[i][j] = maxfn(0, (Nnot - 1) * Sil2[i][j] * (eta[j]) / epsilon[j] - epsilon[j] * num_obs[i][j] / 2);   //R도 업데이트
                    //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];    //평행선 continuation
                }

            }


            // check whether the decision is correct
            int cd_for_one_threshold = 1;   // 지금 system for문 내부임
            for (int d = 0; d < NumThreshold; d++) {    //threshold에 대해서 검사
                for (int j = 0; j < NumConstraint; j++) {   //constraint에 대해서 검사
                    //mean_value[i][j] = sumY[j] / total_obs;
                    if (mean_value[i][j] <= q[d][j] - epsilon[j]) {     //tolerance level을 뺀 q보다 mean value가 작다면,
                        if (Z[i][j][d] == 1) {      // feasible하게 판단했다면,
                            cd_for_one_threshold *= 1;  // correct decision
                        }
                        else {
                            cd_for_one_threshold *= 0;  // incorrect decision
                        }
                    }
                    else if (mean_value[i][j] >= q[d][j] + epsilon[j]) {    //tolerance level을 더한 q보다 mean value가 크다면,
                        if (Z[i][j][d] == 0) {  //만약 infeasible하게 판단했다면, 
                            cd_for_one_threshold *= 1;  //correct decision
                        }
                        else {
                            cd_for_one_threshold *= 0;  //incorrect
                        }
                    }
                    else {
                        cd_for_one_threshold *= 1;  //Acceptable하니깐 그냥 correct decision이라고 하자.
                    }
                }
            }

            if (cd_for_one_threshold == 1) {    //만약 correct decision이면,
                final_cd *= 1;  //최종 correct decision
            }
            else {
                final_cd *= 0;  //최종 incorrect
            }

        }

        if (final_cd == 1) {    //iteration에서 correct decision의 개수를 세주자
            correct_decision++;
        }

        overall_obs += total_obs;   //각 시스템마다 했던 obs 갯수를 다 더함

    }

    printf("Overall: %.10f\n", correct_decision / NumMacro);
    printf("Overall: %.4f\n", overall_obs / NumMacro);

    return 0;
}

double MRG32k3a() //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11 - a13n * s10;
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10 = s11;   s11 = s12;   s12 = p1;

    // Component 2
    p2 = a21 * s22 - a23n * s20;
    k = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20 = s21;   s21 = s22;   s22 = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm) + 0.000001;
}

double generate_Bernoulli(int numConstraint) {

    int numBatches = 1;
    double p = 0.85;

    for (int i = 0; i < NumSys; i++) {
        for (int j = 0; j < numConstraint; j++) {
            int successes = 0;
            for (int k = 0; k < numBatches; k++) {
                if (MRG32k3a() <= p) {
                    successes++;
                }
                observations[i][j] = static_cast<double>(successes) / numBatches;
            }
        }

        return 0;
    }
}

double configuration(void) {


    for (int i = 0; i < NumSys; i++) {
        for (int j = 0; j < NumConstraint; j++) {
            mean_value[i][j] = 0.85;
            ON[i][j] = 1;
            for (int d = 0; d < NumThreshold; d++) {
                ON_l[i][j][d] = 1;
            }
        }
    }
    // Single system
    // system_info[0] = 1;
    for (int j = 0; j < NumConstraint; j++) {
        epsilon[j] = 0.01;
    }
    
    //one threshold
    q[0][0] = 0.86;

    //two threshold
    /*q[0][0] = 0.8;
    q[1][0] = 0.9;*/

    //20 thresholds
    /*for (int j = 0; j < NumConstraint; j++) {
        epsilon[j] = 0.01;
    }

    for (int d = 0; d < NumThreshold; d++) {
        if (d <= 9) {
            q[d][0] = 0.75 + d * epsilon[0];
        }
        else {
            q[d][0] = 0.85 + (d - 9) * epsilon[0];
        }
    }*/

    return 0;
}

double maxfn(double x, double y)
{
    if (x > y) return x;
    else return y;
}

double minfn(double x, double y)
{
    if (x < y) return x;
    else return y;
}
