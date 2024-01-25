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
#include <random>
#include <cmath>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define Nnot	20
#define NumMacro 1000
#define NumSys	1
#define NumConstraint	2
#define NumThreshold	10
#define probability1 0.85
#define probability2 0.6
#define NumBatch    100
#define Theta   1.5

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

double normal(double rmean, double rvar);
double configuration(void);
double generate_Bernoulli(int numConstraint, int case_index);
int read_chol_matrix(void);

double mean_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];
double q[NumThreshold][NumConstraint];
double qU[NumThreshold][NumConstraint];
double qL[NumThreshold][NumConstraint];
double newq[NumThreshold][NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];


int main()
{
    // modify read_chol_matrix() depends on the number of constranits and the correlation between constraints
    read_chol_matrix();

    double total_obs;   //총 obs 갯수
    double final_cd;    //correct decision 인지 아닌지
    double correct_decision;    // 최종 correct decision 갯수
    double overall_obs;     // 최종 obs 갯수

    correct_decision = 0;
    overall_obs = 0;

    // 1 threshold. 95%
    //double alpha = 0.05;
    // >= 2 threshold
    double alpha = 0.025;

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

    // Section 6.1 single system with single constraint
    //eta[0] = 0.1854;

    // Section 6.2 single system with multiple constraints
    /*for (int j = 0; j < NumConstraint; j++) {
        eta[j] = 0.3119;
    }*/

    // Section 6.3 multiple systems with multiple constraints
    /*for (int j=0; j<NumConstraint; j++) {
        eta[j] = 0.5520;
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
                Sil2[i][j] = 0;     //Nnot까지 i 시스템의 분산
                num_obs[i][j] = 0;      //i번째 시스템의 j 번째 제약에서 갯수
            }
        }

        for (int i = 0; i < NumSys; i++) {

            // generate initial samples
            double sumY[NumConstraint];     //obs 값들의 합
            double sum_squareY[NumConstraint];      //obs값들의 sum square
            //double sumYq[NumConstraint];    //obs 값 - q의 값들의 합
            //double avgY[NumConstraint];    // Nnot 까지 obs 값들의 평균

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] = 0;
                sum_squareY[j] = 0;
                //    sumYq[j] = 0;
                //    avgY[j] = 0;
            }

            for (int n = 0; n < Nnot; n++) {    //Nnot 만큼 최초 obs 생성
                generate_Bernoulli(NumConstraint, system_info[i]);  // obs 하나 만듦 (베르누이 분포에서 평균)
                total_obs += 1;     // constraint와 상관없이 시스템의 총 obs 갯수 모든 제약조건에 대해 같이 증가함.


                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    sum_squareY[j] += observations[i][j] * observations[i][j];
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                }



                for (int j = 0; j < NumConstraint; j++) {
                    surviveThreshold[j] = NumThreshold;     //threshold 개수로 survive한거 정의 하고 시작
                    //printf("Overall: %.10f\n", sumY[j] / num_obs[i][j]);
                }
            }

            // find continuation region. 일단 Nnot까지 한거에 대해서 구함
            for (int j = 0; j < NumConstraint; j++) {
                //avgY[j] = sumY[j] / Nnot;
                //Sil2[i][j] = (sum_squareY[j] - (sumY[j] * sumY[j]) / Nnot) / (Nnot - 1);
                Sil2[i][j] = (sum_squareY[j] / (Nnot - 1)) -(sumY[j] / (Nnot - 1)) * (sumY[j] / Nnot);
                R[i][j] = maxfn(0, (Nnot - 1) * Sil2[i][j] * (eta[j]) / epsilon[j] - epsilon[j] * num_obs[i][j] / 2);
                //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];    //평행선 continuation
            }

            while (surviveConstraint != 0) {

                for (int j = 0; j < NumConstraint; j++) {   //모든 const에 대해서 아래 검사

                    if (ON[i][j] == 1) {    //만약 그 시스템에서 그 constraint가 아직 검사가 안됐다면

                        for (int d = 0; d < NumThreshold; d++) {    //threshold까지 내려가서 검사
                            newq[d][j] = 0;
                            qU[d][j] = q[d][j]*Theta/(q[d][j]*Theta + (1-q[d][j]));
                            qL[d][j] = q[d][j]/(q[d][j] + (1-q[d][j])*Theta);
                            newq[d][j] = (qU[d][j] + qL[d][j])/2;
                            if (ON_l[i][j][d] == 1) {   //검사 안됐었으면

                                if ((sumY[j] + R[i][j]) / num_obs[i][j] <= newq[d][j]) {     // feasible 조건
                                    Z[i][j][d] = 1;     // i 시스템의 j 제약의 d번째 threshold에 대해 feasible
                                    ON_l[i][j][d] = 0;  // 해당 threshold를 검사한거로 변경
                                    surviveThreshold[j] -= 1;
                                }

                                else if ((sumY[j] - R[i][j]) / num_obs[i][j] >= newq[d][j]) {   // infeasible 조건
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

                generate_Bernoulli(NumConstraint, system_info[i]);      //만약 위 조건을 만족 안했었다면 obs 하나 더 생성
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
                    newq[d][j] = 0;
                    qU[d][j] = q[d][j]*Theta/(q[d][j]*Theta + (1-q[d][j]));
                    qL[d][j] = q[d][j]/(q[d][j] + (1-q[d][j])*Theta);
                    newq[d][j] = (qU[d][j] + qL[d][j])/2;
                    //mean_value[i][j] = sumY[j] / total_obs;
                    if (mean_value[i][j] <= newq[d][j] - epsilon[j]) {     //tolerance level을 뺀 q보다 mean value가 작다면,
                        if (Z[i][j][d] == 1) {      // feasible하게 판단했다면,
                            cd_for_one_threshold *= 1;  // correct decision
                        }
                        else {
                            cd_for_one_threshold *= 0;  // incorrect decision
                        }
                    }
                    else if (mean_value[i][j] >= newq[d][j] + epsilon[j]) {    //tolerance level을 더한 q보다 mean value가 크다면,
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
    printf("Overall: %.4f\n", (overall_obs * NumBatch) / NumMacro);

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

double normal(double rmean, double rvar)
/* return normal random variable with mean rmean and variance rvar. */
// this is modified for Fully Sequential Procedure with CRN
{
    double V1 = 0, V2 = 0, W = 2, Y = 0, X1 = 0;
    do {
        V1 = 2 * MRG32k3a() - 1;
        V2 = 2 * MRG32k3a() - 1;
        W = pow(V1, 2) + pow(V2, 2);

    } while (W > 1);

    Y = sqrt((-2.00 * log(W)) / W);
    X1 = rmean + sqrt(rvar) * V1 * Y;
    return X1;
}

int read_chol_matrix() {
    double c1, c2, c3, c4, c5;
    char ch;

    // change the input file depending on the correlation between constraints
    std::ifstream myfile("cholMatrix_rho-0.5.txt");
    if (myfile.is_open()) {
        int case_counter = 0;
        int pair_counter = 0;
        while (myfile >> c1 >> ch >> c2 >> ch >> c3 >> ch >> c4 >> ch >> c5)
        {

            // modify the following part depending on the number of constriants considered
            // the current code can handle up to five constraints
            chol_matrix[case_counter][pair_counter][0] = c1;
            chol_matrix[case_counter][pair_counter][1] = c2;
            //chol_matrix[case_counter][pair_counter][2] = c3;
            //chol_matrix[case_counter][pair_counter][3] = c4;
            //chol_matrix[case_counter][pair_counter][4] = c5;
            pair_counter++;
            if (pair_counter == NumConstraint) {
                case_counter++;
                pair_counter = 0;
            }
        }
        myfile.close();
    }
    return 0;
}

//NORTA로 correlated observations 생성
double generate_Bernoulli(int numConstraint, int case_index) {

    std::vector<int> numBatches(numConstraint);
    std::vector<double> p(numConstraint);
    std::vector<double> vari(numConstraint);

    numBatches[0] = NumBatch;
    //numBatches[1] = 10;
    p[0] = probability1;
    p[1] = probability2;
    //vari[0] = p[0] * (1 - p[0]) / numBatches[0];
    //vari[1] = p[1] * (1 - p[1]) / numBatches[0];


    // C matrix is directly input from read_chol_matrix()
    // choose case_index to denote CV, IV or DV
    // observation = mu + CZ
    std::vector<std::vector<double>> C(numConstraint, std::vector<double>(numConstraint));

    switch (case_index) {
    case 1:  // constant variance (CV)

        for (int i = 0; i < numConstraint; i++) {
            for (int j = 0; j < numConstraint; j++) {
                C[i][j] = chol_matrix[0][i][j];

            }
        }
        break;

    case 2:  // increasing variance (IV)

        for (int i = 0; i < numConstraint; i++) {
            for (int j = 0; j < numConstraint; j++) {
                C[i][j] = chol_matrix[1][i][j];
            }
        }
        break;

    case 3:  // decreasing variance (DV)

        for (int i = 0; i < numConstraint; i++) {
            for (int j = 0; j < numConstraint; j++) {
                C[i][j] = chol_matrix[2][i][j];
            }
        }
        break;
    }

    for (int i = 0; i < NumSys; i++) {
        std::vector<double> std_normal(numConstraint);
        std::vector<double> correlated_normal(numConstraint);
        std::vector<int> successes(numConstraint);
        std::vector<double> x_value(numConstraint);

        successes[0] = 0;
        successes[1] = 0;
        boost::math::normal_distribution<> standard_normal;
        x_value[0] = boost::math::quantile(standard_normal, p[0]);
        x_value[1] = boost::math::quantile(standard_normal, p[1]);

        for (int k = 0; k < NumBatch; k++) {
            // Generate independent standard normal random variables
            for (int l = 0; l < numConstraint; l++) {
                std_normal[l] = normal(0, 1);
            }

            // Calculate correlated normal random variables and convert to Bernoulli
            for (int j = 0; j < numConstraint; j++) {
                correlated_normal[j] = 0;
                for (int m = 0; m < numConstraint; m++) {
                    correlated_normal[j] += C[j][m] * std_normal[m];
                }

                if (correlated_normal[j] < x_value[j]) {
                    successes[j]++;
                }
            }
        }

        // Divide the number of successes by the number of batches for each constraint
        for (int j = 0; j < numConstraint; j++) {
            observations[i][j] = static_cast<double>(successes[j]) / NumBatch;
            //printf("obs: % .2f\n", observations[i][j]);
        }
    }

    return 0;
}

double configuration(void) {

    for (int i = 0; i < NumSys; i++) {
        for (int j = 0; j < NumConstraint; j++) {
            mean_value[i][0] = probability1;
            mean_value[i][1] = probability2;
            ON[i][j] = 1;
            for (int d = 0; d < NumThreshold; d++) {
                ON_l[i][j][d] = 1;
            }
        }
    }
    // Single system
    system_info[0] = 1;
    // for (int j = 0; j < NumConstraint; j++) {
    //     epsilon[0] = 0.0205;
    //     epsilon[1] = 0.0402;
    // }


    double maxthreshold[NumConstraint];
    maxthreshold[0] = 0.95;
    maxthreshold[1] = 0.95;

    // two constraints and two thresholds
    // q[0][0] = 0.84;
    // q[1][0] = 0.86;
    // q[0][1] = 0.55;
    // q[1][1] = 0.65;

    // 10 constraints and two thresholds
    q[0][0] = 0.5;
    q[1][0] = 0.55;
    q[2][0] = 0.6;
    q[3][0] = 0.65;
    q[4][0] = 0.7;
    q[5][0] = 0.75;
    q[6][0] = 0.8;
    q[7][0] = 0.85;
    q[8][0] = 0.9;
    q[9][0] = 0.95;
    q[0][1] = 0.5;
    q[1][1] = 0.55;
    q[2][1] = 0.6;
    q[3][1] = 0.65;
    q[4][1] = 0.7;
    q[5][1] = 0.75;
    q[6][1] = 0.8;
    q[7][1] = 0.85;
    q[8][1] = 0.9;
    q[9][1] = 0.95;

    double theta[NumConstraint];
    double UB, LB;

    theta[0] = Theta;
    theta[1] = Theta;

    for (int j = 0; j < NumConstraint; j++) {
        UB = maxthreshold[j]*theta[j]/(maxthreshold[j]*theta[j] + (1-maxthreshold[j]));
        LB = maxthreshold[j]/(maxthreshold[j] + (1-maxthreshold[j])*theta[j]);
        epsilon[j] = (UB-LB)/2;
        //printf("epsilon: %.4f\n", epsilon[j]);
    }

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
