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
    double total_obs;
    double final_cd;
    double correct_decision;
    double overall_obs;

    correct_decision = 0;
    overall_obs = 0;

    double eta[NumConstraint];

    eta[0] = 0.1854;
    //eta[0] = 0.1854;


    for (int l = 0; l < NumMacro; l++) {

        configuration();
        total_obs = 0;
        final_cd = 1;

        double num_obs[NumSys][NumConstraint];
        double R[NumSys][NumConstraint];
        double Sil2[NumSys][NumConstraint];

        for (int i = 0; i < NumSys; i++) {
            for (int j = 0; j < NumConstraint; j++) {
                Sil2[i][j] = 0;
                num_obs[i][j] = 0;
            }
        }

        for (int i = 0; i < NumSys; i++) {

            // generate initial samples
            double sumY[NumConstraint];
            double sum_squareY[NumConstraint];
            double sumYq[NumConstraint];

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] = 0;
                sum_squareY[j] = 0;
                sumYq[j] = 0;
            }

            for (int n = 0; n < Nnot; n++) {
                generate_Bernoulli(NumConstraint);
                total_obs += 1;

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    sum_squareY[j] += observations[i][j] * observations[i][j];
                    sumYq[j] += observations[i][j] - q[0][0];
                    num_obs[i][j] += 1;
                }

                for (int j = 0; j < NumConstraint; j++) surviveThreshold[j] = NumThreshold;

            }

            // find continuation region
            for (int j = 0; j < NumConstraint; j++) {
                Sil2[i][j] = (sum_squareY[j] / (Nnot - 1)) - (sumY[j] / (Nnot - 1)) * (sumY[j] / Nnot);
                R[i][j] = maxfn(0, (Nnot - 1) * Sil2[i][j] * (eta[j]) / epsilon[j] - epsilon[j] * num_obs[i][j] / 2);
                //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];
            }

            while (surviveConstraint != 0) {

                for (int j = 0; j < NumConstraint; j++) {

                    if (ON[i][j] == 1) {

                        for (int d = 0; d < NumThreshold; d++) {

                            if (ON_l[i][j][d] == 1) {

                                if ((sumY[j] + R[i][j]) / num_obs[i][j] <= q[d][j]) {
                                    Z[i][j][d] = 1;
                                    ON_l[i][j][d] = 0;
                                    surviveThreshold[j] -= 1;
                                }

                                if ((sumY[j] - R[i][j]) / num_obs[i][j] >= q[d][j]) {
                                    Z[i][j][d] = 0;
                                    ON_l[i][j][d] = 0;
                                    surviveThreshold[j] -= 1;
                                }
                            }

                        }

                        if (surviveThreshold[j] == 0) {
                            ON[i][j] = 0;
                            surviveConstraint -= 1;
                        }

                    }

                }

                if (surviveConstraint == 0) break;

                generate_Bernoulli(NumConstraint);
                total_obs += 1;

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    num_obs[i][j] += 1;
                    R[i][j] = maxfn(0, (Nnot - 1) * Sil2[i][j] * (eta[j]) / epsilon[j] - epsilon[j] * num_obs[i][j] / 2);
                    //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];
                }

            }

            // check whether the decision is correct
            int cd_for_one_threshold = 1;
            for (int d = 0; d < NumThreshold; d++) {
                for (int j = 0; j < NumConstraint; j++) {
                    if (mean_value[i][j] <= q[d][j] - epsilon[j]) {
                        if (Z[i][j][d] == 1) {
                            cd_for_one_threshold *= 1;
                        }
                        else {
                            cd_for_one_threshold *= 0;
                        }
                    }
                    else if (mean_value[i][j] >= q[d][j] + epsilon[j]) {
                        if (Z[i][j][d] == 0) {
                            cd_for_one_threshold *= 1;
                        }
                        else {
                            cd_for_one_threshold *= 0;
                        }
                    }
                    else {
                        cd_for_one_threshold *= 1;
                    }
                }
            }

            if (cd_for_one_threshold == 1) {
                final_cd *= 1;
            }
            else {
                final_cd *= 0;
            }

        }

        if (final_cd == 1) {
            correct_decision++;
        }

        overall_obs += total_obs;

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

        q[0][0] = 0.84;
        q[1][0] = 0.86;

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