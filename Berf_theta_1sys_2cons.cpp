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
#include <cmath>
#include <random>
#include <list>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
using namespace std;

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)

#define NumMacro 1000
#define NumSys	1
#define NumConstraint	2
#define NumThreshold	2

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

double normal(double rmean, double rvar);
double configuration(void);
//double generate_multiNormal(int numConstraint, int case_index);
double generate_Bernoulli(int numConstraint, int case_index);
int read_chol_matrix(void);
float correlationCoefficient(std::vector<double> X, std::vector<double> Y, int n);

double mean_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];
double H[NumSys][NumConstraint];
double dummies[NumThreshold][NumConstraint];
double q[NumThreshold][NumConstraint];
double qL[NumThreshold][NumConstraint];
double qU[NumThreshold][NumConstraint];
double alpha;
double theta[NumConstraint];
double thetas[NumThreshold][NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];

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

    std::ifstream myfile("cholMatrix_rho0.txt");
    if (myfile.is_open()) {
        int case_counter = 0;
        int row_counter = 0;  // renamed to make the purpose more clear
        int max_rows = 2;    // we're only interested in the first 2 rows
        while (myfile >> c1 >> ch >> c2 >> ch >> c3 >> ch >> c4 >> ch >> c5)
        {
            // We're only interested in the first 2*2 matrix in each 5*5 matrix
            if (row_counter < max_rows) {
                chol_matrix[case_counter][row_counter][0] = c1;
                chol_matrix[case_counter][row_counter][1] = c2;
                //chol_matrix[case_counter][row_counter][2] = c3;
                // you might need to initialize the rest of the elements in this row to some default value
            }
            row_counter++;
            if (row_counter == 5) {  // we've read a complete 5*5 matrix
                case_counter++;
                row_counter = 0;
            }
            if (case_counter == 3) {  // we've read all 3 5*5 matrices
                break;
            }
        }
        myfile.close();
    }
    return 0;
}

double generate_Bernoulli(int numConstraint, int case_index) {

    std::vector<int> numBatches(numConstraint);
    std::vector<double> truep(numConstraint);
    //std::vector<double> vari(numConstraint);

    numBatches[0] = 1;
    //numBatches[1] = 10;
    truep[0] = 0.1;
    truep[1] = 0.4;
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
                //printf("C: % .2f\n", C[i][j]);
            }
        }
        break;

    case 2:  // increasing variance (IV)

        for (int i = 0; i < numConstraint; i++) {
            for (int j = 0; j < numConstraint; j++) {
                C[i][j] = chol_matrix[1][i][j];
                //printf("C: % .2f\n", C[i][j]);
            }
        }
        break;

    case 3:  // decreasing variance (DV)

        for (int i = 0; i < numConstraint; i++) {
            for (int j = 0; j < numConstraint; j++) {
                C[i][j] = chol_matrix[2][i][j];
                //printf("C: % .2f\n", C[i][j]);
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
        x_value[0] = boost::math::quantile(standard_normal, truep[0]);
        x_value[1] = boost::math::quantile(standard_normal, truep[1]);

        for (int k = 0; k < numBatches[0]; k++) {
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

                for (int d = 0; d < NumThreshold; d++) {
                    double rn = MRG32k3a();
                    dummies[d][j] = 0;
                    if (rn <= q[d][j]) {
                        dummies[d][j] = 1;
                    }
                }
            }
        }

        // Divide the number of successes by the number of batches for each constraint
        for (int j = 0; j < numConstraint; j++) {
            observations[i][j] = static_cast<double>(successes[j]) / numBatches[0];
            //printf("obs: % .2f\n", observations[i][j]);
        }
    }

    return 0;
}

double configuration(void) {

    for (int i = 0; i < NumSys; i++) {
        for (int j = 0; j < NumConstraint; j++) {
            mean_value[i][0] = 0.1;
            mean_value[i][1] = 0.4;
            ON[i][j] = 1;
            for (int d = 0; d < NumThreshold; d++) {
                ON_l[i][j][d] = 1;
            }
        }
    }
    // Single system
    system_info[0] = 1;

    /*for (int j = 0; j < NumConstraint; j++) {
        epsilon[0] = 0.05;
        epsilon[1] = 0.05;
    }*/

    // two constraints and two thresholds
    q[0][0] = 0.09;
    q[1][0] = 0.11;
    //q[2][0] = 0.70;
    q[0][1] = 0.35;
    q[1][1] = 0.45;
    //q[2][1] = 0.30;

    for (int j = 0; j < NumConstraint; j++) {
        theta[j] = 1.5;
    }

    /*for (int j = 0; j < NumConstraint; j++) {
        double arr[NumThreshold];
        for (int d = 0; d < NumThreshold; d++) {
            qL[d][j] = q[d][j] - epsilon[j];
            qU[d][j] = q[d][j] + epsilon[j];
            thetas[d][j] = std::min(qU[d][j] / (1.0 - qU[d][j]) * (1.0 - q[d][j]) / q[d][j], q[d][j] / (1 - q[d][j]) * (1 - qL[d][j]) / qL[d][j]);
            arr[d] = thetas[d][j];
        }
        theta[j] = *std::min_element(std::begin(arr), std::end(arr));
    }*/

    return 0;
}

float correlationCoefficient(std::vector<double> X, std::vector<double> Y, int n)
{

    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }

    // use formula for calculating correlation coefficient.
    float corr = (float)(n * sum_XY - sum_X * sum_Y)
        / sqrt((n * squareSum_X - sum_X * sum_X)
            * (n * squareSum_Y - sum_Y * sum_Y));

    return corr;
}

struct f {
    //int idx;
    //f(int index) : idx(index) {}
    // 1 threshold. 95%
    double alpha = 0.05;
    // >= 2 threshold
    //double alpha = 0.025;
    double k = NumSys;
    double s = NumConstraint;
    double beta = alpha / (2 * k * s);

    double operator()(double ell) {
        return (1 / (1 + std::pow(theta[0], ell)) - beta);
    }
};

double fsolve() {
    const boost::uintmax_t maxit = 1000;
    boost::uintmax_t it = maxit;
    const double guess = 1.0;
    const double min = 0.0;
    const double max = 2.0;
    boost::math::tools::eps_tolerance<double> tol((std::numeric_limits<double>::digits * 3) / 2);
    std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(f(), guess, 2.0, false, tol, it);
    return r.first + (r.second - r.first) / 2;
};

int main()
{
    read_chol_matrix();

    double total_obs;   //총 obs 갯수 in an iter
    double final_cd;    //correct decision 인지 아닌지
    double correct_decision;    // 최종 correct decision 갯수
    double overall_obs;     // 최종 obs 갯수
    double overall_corr;

    correct_decision = 0;
    overall_obs = 0;
    overall_corr = 0;

    for (int l = 0; l < NumMacro; l++) {    // iteration 시작. 

        configuration();    //최초설정
        total_obs = 0;      //우선 obs 갯수 0으로 시작
        final_cd = 1;       //일단 cd를 1로 시작
        std::vector<double> X;
        std::vector<double> Y;

        double num_obs[NumSys][NumConstraint];

        for (int i = 0; i < NumSys; i++) {
            for (int j = 0; j < NumConstraint; j++) {
                H[i][j] = std::ceil(fsolve());
                //printf("theta: %.10f\n", theta[j]);
                //printf("H: %.10f\n", H[i][j]);
            }
        }

        for (int i = 0; i < NumSys; i++) {

            // generate initial samples
            double sumY[NumConstraint];     //obs 값들의 합
            double sumI[NumThreshold][NumConstraint];     //obs 값들의 합

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] = 0;
                for (int d = 0; d < NumThreshold; d++) {
                    sumI[d][j] = 0;
                }
            }

            //최초 1회 시행
            generate_Bernoulli(NumConstraint, system_info[i]);  // obs 하나 만듦 (베르누이 분포에서 평균)
            total_obs += 1;     // constraint와 상관없이 시스템의 총 obs 갯수 모든 const에 대해 같이 증가함.

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] += observations[i][j];
                for (int d = 0; d < NumThreshold; d++) {
                    sumI[d][j] += dummies[d][j];
                }
                num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                if (j == 0) {
                    X.push_back(observations[i][j]);
                }
                else if (j == 1) {
                    Y.push_back(observations[i][j]);
                }
            }

            for (int j = 0; j < NumConstraint; j++) surviveThreshold[j] = NumThreshold;     //threshold 개수로 survive한거 정의 하고 시작

            while (surviveConstraint != 0) {

                for (int j = 0; j < NumConstraint; j++) {   //모든 const에 대해서 아래 검사

                    if (ON[i][j] == 1) {    //만약 그 시스템에서 그 constraint가 아직 검사가 안됐다면

                        for (int d = 0; d < NumThreshold; d++) {    //threshold까지 내려가서 검사

                            if (ON_l[i][j][d] == 1) {   //검사 안됐었으면

                                if (sumY[j] - sumI[d][j] <= -H[i][j]) {     // 첫번째 조건
                                    Z[i][j][d] = 1;     //
                                    ON_l[i][j][d] = 0;  // 해당 threshold를 검사한거로 변경
                                    surviveThreshold[j] -= 1;
                                }

                                else if (sumY[j] - sumI[d][j] >= H[i][j]) {   // 두번째 조건
                                    Z[i][j][d] = 0;     //
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
                total_obs += 1;     //total obs도 하나 늘리기

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[d][j] += dummies[d][j];
                    }
                    num_obs[i][j] += 1;
                    if (j == 0) {
                        X.push_back(observations[i][j]);
                    }
                    else if (j == 1) {
                        Y.push_back(observations[i][j]);
                    }
                }

            }

            // check whether the decision is correct
            int cd_for_one_threshold = 1;   // 지금 system for문 내부임
            for (int d = 0; d < NumThreshold; d++) {    //threshold에 대해서 검사
                for (int j = 0; j < NumConstraint; j++) {   //constraint에 대해서 검사
                    //mean_value[i][j] = sumY[j] / total_obs;
                    if (mean_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * theta[j]))) {     //tolerance level을 뺀 q보다 mean value가 작다면,
                        if (Z[i][j][d] == 1) {      // feasible하게 판단했다면,
                            cd_for_one_threshold *= 1;  // correct decision
                        }
                        else {
                            cd_for_one_threshold *= 0;  // incorrect decision
                        }
                    }
                    else if (mean_value[i][j] >= ( q[d][j] * theta[j] / ( (1 - q[d][j]) + q[d][j] * theta[j]) ) ) {    //tolerance level을 더한 q보다 mean value가 크다면,
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

            //int n = X.size();
            //double correlation = correlationCoefficient(X, Y, n);
            //printf("corr: %.5f\n", correlation);
            //overall_corr += correlation;

        }

        if (final_cd == 1) {    //최종 1이면 correct decision의 개수 하나 올려주기
            correct_decision++;
        }

        overall_obs += total_obs;   //각 iter마다 했던 obs 갯수를 다 더함

    }

    printf("Overall PCD: %.10f\n", correct_decision / NumMacro);
    printf("Overall OBS: %.4f\n", overall_obs / NumMacro);
    printf("Overall correlation: %.4f\n", overall_corr / NumMacro);

    return 0;
}
