#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <random>
#include <list>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace std;

#include <random>

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define NumMacro 100
#define NumSys	77
#define NumConstraint	2
#define NumThreshold	4
#define NumPass 2
#define Num_s 20
#define Num_S 20
#define Theta   1.5
#define CRN 1   //If you use CRN, then this is 1. Or this is 0.
#define UR 1   //If you use U_r, then this is 1. Or this is 0.
#define DOMBERF 1   //If you do mberf, then this is 1. Or this is 0 (Only BeRF is implemented.)

// inputs for Generate R(0,1) by L'ecuyer (1997)
#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

// choices of seeds for Generate R(0,1) by L'ecuyer (1997)
//double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;
double  ss10 = 12345, ss11 = 12345, ss12 = 12345, ss20 = 12345, ss21 = 12345, ss22 = 12345;

double s10[NumSys][NumConstraint];
double s11[NumSys][NumConstraint];
double s12[NumSys][NumConstraint];
double s20[NumSys][NumConstraint];
double s21[NumSys][NumConstraint];
double s22[NumSys][NumConstraint];
double s10_y[NumSys][NumConstraint];
double s11_y[NumSys][NumConstraint];
double s12_y[NumSys][NumConstraint];
double s20_y[NumSys][NumConstraint];
double s21_y[NumSys][NumConstraint];
double s22_y[NumSys][NumConstraint];

double s10a[NumSys];
double s11a[NumSys];
double s12a[NumSys];
double s20a[NumSys];
double s21a[NumSys];
double s22a[NumSys];
double s10_ya[NumSys];
double s11_ya[NumSys];
double s12_ya[NumSys];
double s20_ya[NumSys];
double s21_ya[NumSys];
double s22_ya[NumSys];

double MRG32k3a(int sys_index, int constraint_index);  //Generate R(0,1) by L'ecuyer (1997)
double MRG32k3a_y(int sys_index, int constraint_index);
double MRG32k3a_y2(int sys_index);
double MRG32k3a2(int sys_index);
double MRG32k3a3(void);
double minfn(double x1, double x2);
double maxfn(double x, double y);
double normal(double rmean, double rvar, int sys_index, int constraint_index);
double poisson2(double pmean);
double configuration(void);
int read_rand_seeds(void);
// double generate_one_obs(int demand_index, int sys_index, int constraint_index);
double generate_one_obs(int demand_index, int sys_index);
int generate_demand(int sys_index, int constraint_index);
// int generate_demand(void);
int write_up(void);
int read_system_true_value(void);
int determine_true_feasibility(void);

int berf(void);
int mberf1(int pass_index);
int mberf2(int pass_index);

double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];
double H[NumConstraint];
double dummies[NumConstraint][NumThreshold];
double q[NumThreshold][NumConstraint];
double qL[NumConstraint][NumThreshold];
double qU[NumConstraint][NumThreshold];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int BeRF_Z[NumSys][NumConstraint][NumThreshold];
int MBeRF_Z[NumSys][NumConstraint][NumThreshold];
double BeRF_total_obs[NumConstraint];
double MBeRF_total_obs[NumConstraint];
double MBeRF_rep_by_pass[NumPass];
int T_index[NumPass][NumThreshold][NumConstraint];
int system_value[NumSys][2];
// double system_true_value[NumSys][2] = {0};
double system_true_value[NumSys][2];
int true_feasibility[NumSys][NumConstraint][NumThreshold];
double single_obs[NumConstraint];
double demand_list[20000000];
double prn_list[20000000];
double sumY[NumSys][NumConstraint];
double sumI[NumSys][NumConstraint][NumThreshold];
double num_obs[NumSys][NumConstraint];
double v_UB[NumSys][NumConstraint];
double v_LB[NumSys][NumConstraint];
int LAST[NumSys][NumConstraint];

double demand_mean = 25;
double order_cost = 3;
double fixed_order_cost = 32;
double holding_cost = 1;
double penalty_cost = 5;

double total_correct_berf;
double total_correct_mberf;
double total_match_decision;
double berf_total;
double mberf_total;
double berf_per_macro;
double mberf_per_macro[NumPass];
double correct_berf;
double correct_mberf;
double mrf_total_rep_per_macro;

FILE *outfile;

int mberf1(int pass_index) {

    mberf_per_macro[pass_index] = 0;

    for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
			ON[i][j] = 1;
			for (int d=0; d<NumThreshold; d++) {
                if (T_index[pass_index][d][j] == 1) {
                    ON_l[i][j][d] = 1;
                } else ON_l[i][j][d] = 0;
			}
		}
	}

    for (int i=0; i<NumSys; i++) {
   	    for (int j=0; j<NumConstraint; j++) {
            num_obs[i][j] = 0;
            v_UB[i][j] = 1000000000000;
            v_LB[i][j] = -1000000000000;
            LAST[i][j] = -2;
   		}
   	}

    for (int i=0; i<NumSys; i++) {

   	    // generate initial samples
        int surviveConstraint = NumConstraint;
        int surviveThreshold[NumConstraint];

        int demand_index = 0;
        int prn_index = 0;

   		for (int j=0; j<NumConstraint; j++) {
   	        sumY[i][j] = 0;
   			for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] = 0;
            }
            surviveThreshold[j] = 0;
   		}

        generate_one_obs(demand_index, i);

        if (UR == 1){
            double PRN = prn_list[prn_index];
            prn_index += 1;
            for (int j = 0; j < NumConstraint; j++) {
                for (int d = 0; d < NumThreshold; d++) {
                    dummies[j][d] = 0;
                    if (PRN <= q[d][j]) {
                        dummies[j][d] = 1;
                    }
                }
            }
        }

        for (int j = 0; j < NumConstraint; j++) {
            
            sumY[i][j] += single_obs[j];
            for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] += dummies[j][d];
            }
            num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
            MBeRF_total_obs[j] += 1;
        }

        demand_index += 12;
        mberf_total += 1;
        MBeRF_rep_by_pass[pass_index] += 1;
        mberf_per_macro[pass_index] += 1;

        // for same number of thresholds per constriant
        for (int j=0; j<NumConstraint; j++) {
            for (int d=0; d<NumThreshold; d++) surviveThreshold[j] += T_index[pass_index][d][j];
        }

   		while (surviveConstraint != 0) {

   			for (int j=0; j<NumConstraint; j++) {

   				// if ((v_UB[i][j] > v_LB[i][j])) {

                    if ((sumY[i][j]+H[j])/num_obs[i][j] < v_UB[i][j]) LAST[i][j] = 0;   // UB is 0
                    if ((sumY[i][j]-H[j])/num_obs[i][j] > v_LB[i][j]) LAST[i][j] = 1;   // LB is 1

                    v_UB[i][j] = minfn(v_UB[i][j], (sumY[i][j]+H[j])/num_obs[i][j]);
                    v_LB[i][j] = maxfn(v_LB[i][j], (sumY[i][j]-H[j])/num_obs[i][j]);

   					for (int d=0; d<NumThreshold; d++) {

                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                if ((sumY[i][j]+H[j])/num_obs[i][j] <= (sumI[i][j][d])/num_obs[i][j]) {
       								MBeRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                else if ((sumY[i][j]-H[j])/num_obs[i][j] >= (sumI[i][j][d])/num_obs[i][j]) {
                                    MBeRF_Z[i][j][d] = 0;
       								ON_l[i][j][d] = 0;
       				                surviveThreshold[j] -= 1;
       							}
       						}
                        }
                    }

                    if (ON[i][j] == 1) {
                        if (surviveThreshold[j] == 0) {
                            ON[i][j] = 0;
                            surviveConstraint -= 1;
                        }
                    }
   				// }
   			}

   			if (surviveConstraint == 0) {
                // printf("MBeRF1 finished\n");
                break;
            }

            generate_one_obs(demand_index, i);

            if (UR == 1){
                double PRN = prn_list[prn_index];
                prn_index += 1;
                for (int j = 0; j < NumConstraint; j++) {
                    for (int d = 0; d < NumThreshold; d++) {
                        dummies[j][d] = 0;
                        if (PRN <= q[d][j]) {
                            dummies[j][d] = 1;
                        }
                    }
                }
            }

   			for (int j=0; j<NumConstraint; j++) {

                // if (v_UB[i][j] > v_LB[i][j]) {
                if (ON[i][j] == 1) {
                    // generate_one_obs(demand_index, i, j);
                    // printf("2nd obs: %.1f\n",single_obs[1]);
                    sumY[i][j] += single_obs[j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                        // printf("2nd dummy obs: %.1f\n",dummies[1][d]);
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    MBeRF_total_obs[j] += 1;    
                }            
                // }
   			}

            demand_index += 12;
            mberf_total += 1;
            MBeRF_rep_by_pass[pass_index] += 1;
            mberf_per_macro[pass_index] += 1;
   		}
   	}

    return 0;
}

int mberf2(int pass_index) {

    mberf_per_macro[pass_index] = 0;

    for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
            ON[i][j] = 0;
            for (int d=0; d<NumThreshold; d++) {
                if (T_index[pass_index][d][j] == 1) {
                    ON[i][j] = 1;
                    ON_l[i][j][d] = 1;
                } else ON_l[i][j][d] = 0;
            }
        }
	}

    for (int i=0; i<NumSys; i++) {
        int demand_index = 10000000;
        int prn_index = 10000000;

        int surviveConstraint = 0;
        int surviveThreshold[NumConstraint];

        for (int j=0; j<NumConstraint; j++) {
            surviveThreshold[j] = 0;
            for (int d=0; d<NumThreshold; d++) {
                if (T_index[pass_index][d][j] == 1) {
                    surviveThreshold[j] += 1;
                }
            }
            if (surviveThreshold[j] > 0) surviveConstraint += 1;
        }

        while (surviveConstraint != 0) {

            for (int j=0; j<NumConstraint; j++) {

                if (ON[i][j] == 1) {

                    for (int d=0; d<NumThreshold; d++) {

                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                // printf("v_LB: %.5f\n", v_LB[i][j]);
                                // printf("v_UB: %.5f\n", v_UB[i][j]);
                                // printf("E(I): %.5f\n", (sumI[i][j][d]/num_obs[i][j]));    

                                if (((v_UB[i][j] <= q[d][j]) & (q[d][j] < v_LB[i][j])) |
                                    ((v_UB[i][j] < q[d][j]) & (q[d][j] <= v_LB[i][j]))) {

                                        if (LAST[i][j] == 0) {
                                            MBeRF_Z[i][j][d] = 0;
                                            ON_l[i][j][d] = 0;
                                            surviveThreshold[j] -= 1;
                                        }

                                        if (LAST[i][j] == 1) {
                                            MBeRF_Z[i][j][d] = 1;
                                            ON_l[i][j][d] = 0;
                                            surviveThreshold[j] -= 1;
                                        }
                                } else if ((v_UB[i][j] <= q[d][j]) & (v_LB[i][j] <= q[d][j])) {
                                    MBeRF_Z[i][j][d] = 1;
                                    ON_l[i][j][d] = 0;
                                    surviveThreshold[j] -= 1;
                                } else if ((v_LB[i][j] >= q[d][j]) & (v_UB[i][j] >= q[d][j])) {
                            //        if (j == 0) printf("%.5f\t%.5f\t%d\n", v_LB[i][j], v_UB[i][j], d);
                                    MBeRF_Z[i][j][d] = 0;
                                    ON_l[i][j][d] = 0;
                                    surviveThreshold[j] -= 1;
                                }
                            }
                        }
                    }

                    //if (j == 0) printf("%.5f\t%.5f\t%d\t%d\n", v_LB[i][j], v_UB[i][j], surviveThreshold[j], surviveConstraint);
                    if (surviveThreshold[j] == 0) {
                        ON[i][j] = 0;
                        surviveConstraint -= 1;
                    }
                }
            }

    //        printf("%d\n", surviveConstraint);
            if (surviveConstraint == 0) {
                // printf("MBeRF2-1 finished\n");
                break;}

            generate_one_obs(demand_index, i);

            if (UR == 1){
                double PRN = prn_list[prn_index];
                prn_index += 1;
                for (int j = 0; j < NumConstraint; j++) {
                    for (int d = 0; d < NumThreshold; d++) {
                        dummies[j][d] = 0;
                        if (PRN <= q[d][j]) {
                            dummies[j][d] = 1;
                        }
                    }
                }
            }

            for (int j=0; j<NumConstraint; j++) {

                if (v_UB[i][j] > v_LB[i][j]) {
                    // generate_one_obs(demand_index, i, j);
                    sumY[i][j] += single_obs[j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    MBeRF_total_obs[j] += 1;  
                }
   			}
            demand_index += 12;
            mberf_total += 1;
            MBeRF_rep_by_pass[pass_index] += 1;
            mberf_per_macro[pass_index] += 1;

            for (int j=0; j<NumConstraint; j++) {

    //            if (j == 0) printf("%.5f\t%.5f\n", v_LB[i][j], v_UB[i][j]);

   				if ((ON[i][j] == 1) & (v_UB[i][j] > v_LB[i][j])) {

                    if ((sumY[i][j]+H[j])/num_obs[i][j] < v_UB[i][j]) LAST[i][j] = 0;   // UB is 0
                    if ((sumY[i][j]-H[j])/num_obs[i][j] > v_LB[i][j]) LAST[i][j] = 1;   // LB is 1

                    v_UB[i][j] = minfn(v_UB[i][j], (sumY[i][j]+H[j])/num_obs[i][j]);
                    v_LB[i][j] = maxfn(v_LB[i][j], (sumY[i][j]-H[j])/num_obs[i][j]);

   					for (int d=0; d<NumThreshold; d++) {

                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                if (v_UB[i][j] <= (sumI[i][j][d]/num_obs[i][j])) {
       								MBeRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                if (v_LB[i][j] >= (sumI[i][j][d]/num_obs[i][j])) {
                                    MBeRF_Z[i][j][d] = 0;
       								ON_l[i][j][d] = 0;
       				                surviveThreshold[j] -= 1;
       							}
       						}
                        }
                    }

                    if (surviveThreshold[j] == 0) {
                        ON[i][j] = 0;
                        surviveConstraint -= 1;
                    }

   				}
        //    printf("%.5f\t%.5f\t%d\n", v_LB[i][j], v_UB[i][j], surviveThreshold[j]);
   			}

   			if (surviveConstraint == 0) {
                // printf("MBeRF2-2 finished\n");
                break;}
        }
    }

    return 0;
}

int berf(void) {

    berf_per_macro = 0;

    for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
			ON[i][j] = 1;
			for (int d=0; d<NumThreshold; d++) {
				ON_l[i][j][d] = 1;
			}
		}
	}

    for (int i=0; i<NumSys; i++) {
   	    for (int j=0; j<NumConstraint; j++) {
            num_obs[i][j] = 0;
   		}
   	}

   	for (int i=0; i<NumSys; i++) {
        int demand_index = 0;
        int prn_index = 0;

   	    // generate initial samples
        int surviveConstraint = NumConstraint;
        int surviveThreshold[NumConstraint];

        for (int j=0; j<NumConstraint; j++) {
   	        sumY[i][j] = 0;
   			for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] = 0;
            }
            surviveThreshold[j] = 0;
   		}

        generate_one_obs(demand_index, i);

        if (UR == 1){
            double PRN = prn_list[prn_index];
            prn_index += 1;
            for (int j = 0; j < NumConstraint; j++) {
                for (int d = 0; d < NumThreshold; d++) {
                    dummies[j][d] = 0;
                    if (PRN <= q[d][j]) {
                        dummies[j][d] = 1;
                    }
                }
            }
        }

        for (int j = 0; j < NumConstraint; j++) {
            // generate_one_obs(demand_index, i, j);
            sumY[i][j] += single_obs[j];
            for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] += dummies[j][d];
            }
            num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
            BeRF_total_obs[j] += 1;
        }

        demand_index += 12;
        berf_total += 1;
        berf_per_macro += 1;

        // for same number of thresholds per constriant
        for (int j=0; j<NumConstraint; j++) surviveThreshold[j] = NumThreshold;

   		while (surviveConstraint != 0) {

   			for (int j=0; j<NumConstraint; j++) {

   				if (ON[i][j] == 1) {

   					for (int d=0; d<NumThreshold; d++) {

   						if (ON_l[i][j][d] == 1) {

                            if ((sumY[i][j]+H[j])/num_obs[i][j] <= (sumI[i][j][d]/num_obs[i][j])) {
   								BeRF_Z[i][j][d] = 1;
   								ON_l[i][j][d] = 0;
   								surviveThreshold[j] -= 1;
                            }

                            else if ((sumY[i][j]-H[j])/num_obs[i][j] >= (sumI[i][j][d]/num_obs[i][j])) {
                                BeRF_Z[i][j][d] = 0;
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

   			if (surviveConstraint == 0) {
                // printf("BeRF finished\n");
                break;
            }

            generate_one_obs(demand_index, i);

            if (UR == 1){
                double PRN = prn_list[prn_index];
                prn_index += 1;
                for (int j = 0; j < NumConstraint; j++) {
                    for (int d = 0; d < NumThreshold; d++) {
                        dummies[j][d] = 0;
                        if (PRN <= q[d][j]) {
                            dummies[j][d] = 1;
                        }
                    }
                }
            }

   			for (int j=0; j<NumConstraint; j++) {
                if (ON[i][j] == 1) {
                    // generate_one_obs(demand_index, i, j);
                    sumY[i][j] += single_obs[j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    BeRF_total_obs[j] += 1;
                }
   			}
            demand_index += 12;
            berf_total += 1;
            berf_per_macro += 1;
   		}

   	}

	return 0;
}

double MRG32k3a(int sys_index, int constraint_index) //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11[sys_index][constraint_index] - a13n * s10[sys_index][constraint_index];
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10[sys_index][constraint_index] = s11[sys_index][constraint_index];
    s11[sys_index][constraint_index] = s12[sys_index][constraint_index];
    s12[sys_index][constraint_index] = p1;

    // Component 2
    p2 = a21 * s22[sys_index][constraint_index] - a23n * s20[sys_index][constraint_index];
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20[sys_index][constraint_index] = s21[sys_index][constraint_index];
    s21[sys_index][constraint_index] = s22[sys_index][constraint_index];
    s22[sys_index][constraint_index] = p2;

    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double MRG32k3a2(int sys_index) //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11a[sys_index] - a13n * s10a[sys_index];
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10a[sys_index] = s11a[sys_index];
    s11a[sys_index] = s12a[sys_index];
    s12a[sys_index] = p1;

    // Component 2
    p2 = a21 * s22a[sys_index] - a23n * s20a[sys_index];
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20a[sys_index] = s21a[sys_index];
    s21a[sys_index] = s22a[sys_index];
    s22a[sys_index] = p2;

    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double MRG32k3a_y(int sys_index, int constraint_index) //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11_y[sys_index][constraint_index] - a13n * s10_y[sys_index][constraint_index];
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10_y[sys_index][constraint_index] = s11_y[sys_index][constraint_index];   
    s11_y[sys_index][constraint_index] = s12_y[sys_index][constraint_index];
    s12_y[sys_index][constraint_index] = p1;

    // Component 2
    p2 = a21 * s22_y[sys_index][constraint_index] - a23n * s20_y[sys_index][constraint_index];
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20_y[sys_index][constraint_index] = s21_y[sys_index][constraint_index];   
    s21_y[sys_index][constraint_index] = s22_y[sys_index][constraint_index];   
    s22_y[sys_index][constraint_index] = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double MRG32k3a_y2(int sys_index) //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11_ya[sys_index] - a13n * s10_ya[sys_index];
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10_ya[sys_index] = s11_ya[sys_index];   
    s11_ya[sys_index] = s12_ya[sys_index];
    s12_ya[sys_index] = p1;

    // Component 2
    p2 = a21 * s22_ya[sys_index] - a23n * s20_ya[sys_index];
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20_ya[sys_index] = s21_ya[sys_index];   
    s21_ya[sys_index] = s22_ya[sys_index];   
    s22_ya[sys_index] = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double MRG32k3a3() //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;

    // Component 1
    p1 = a12 * ss11 - a13n * ss10;
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    ss10 = ss11;   ss11 = ss12;   ss12 = p1;

    // Component 2
    p2 = a21 * ss22 - a23n * ss20;
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    ss20 = ss21;   ss21 = ss22;   ss22 = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double normal(double rmean, double rvar, int sys_index, int constraint_index)
/* return normal random variable with mean rmean and variance rvar. */
// this is modified for Fully Sequential Procedure with CRN
{
	double V1 = 0, V2 = 0, W = 2, Y = 0, X1 = 0;
	do {
		V1 = 2 * MRG32k3a(sys_index, constraint_index) - 1;
		V2 = 2 * MRG32k3a(sys_index, constraint_index) - 1;
     	W = pow(V1,2) + pow(V2,2);

	} while (W > 1);

	Y = sqrt( (-2.00 * log(W))/W );
	X1 = rmean + sqrt(rvar) * V1 * Y;
	return X1;
}

double normal2(double rmean, double rvar, int sys_index)
/* return normal random variable with mean rmean and variance rvar. */
// this is modified for Fully Sequential Procedure with CRN
{
	double V1 = 0, V2 = 0, W = 2, Y = 0, X1 = 0;
	do {
		V1 = 2 * MRG32k3a2(sys_index) - 1;
		V2 = 2 * MRG32k3a2(sys_index) - 1;
     	W = pow(V1,2) + pow(V2,2);

	} while (W > 1);

	Y = sqrt( (-2.00 * log(W))/W );
	X1 = rmean + sqrt(rvar) * V1 * Y;
	return X1;
}

double poisson(double lam, int sys_index, int constraint_index) {
    double a, b;
    int i;
    a=exp(-lam);
    b=1;
    i=0;

    while(1) {
        b=b*MRG32k3a(sys_index, constraint_index);
        if( b<a) {
            return i;
            break;
        }
        i++;
    }
}

double poisson2(double lam, int sys_index) {
    double a, b;
    int i;
    a=exp(-lam);
    b=1;
    i=0;

    while(1) {
        b=b*MRG32k3a2(sys_index);
        if( b<a) {
            return i;
            break;
        }
        i++;
    }
}

double poisson3(double lam) {
    double a, b;
    int i;
    a=exp(-lam);
    b=1;
    i=0;

    while(1) {
        b=b*MRG32k3a3();
        if( b<a) {
            return i;
            break;
        }
        i++;
    }
}

int read_system_true_value(void) {

  double prob;
  double prob2;

  std::ifstream myfile1 ("real_true_value.txt");
  if (myfile1.is_open()) {
    int system_counter = 0;
    while ( system_counter < NumSys && myfile1 >> prob ) //>> ch >> expected_cost )
    {
      system_true_value[system_counter][0] = prob;
      //system_true_value[system_counter][1] = expected_cost;
      system_counter++;
    }
    myfile1.close();
  }

  std::ifstream myfile2 ("real_true_value2.txt");
  if (myfile2.is_open()) {
    int system_counter = 0;
    while ( system_counter < NumSys && myfile2 >> prob2 ) //>> ch >> expected_cost )
    {
      system_true_value[system_counter][1] = prob2;
      //system_true_value[system_counter][1] = expected_cost;
      system_counter++;
    }
    myfile2.close();
  }

//   for (int i = 0; i < NumSys; i++) printf("sys_true: %.10f\t%.10f\n", system_true_value[i][0], system_true_value[i][1]);
  return 0;
}

int determine_true_feasibility(void) {

  for (int i=0; i<NumSys; i++) {
    for (int j=0; j<NumConstraint; j++) {
      for (int d=0; d<NumThreshold; d++) {
        if (system_true_value[i][j] <= q[d][j]) true_feasibility[i][j][d] = 1;
        else true_feasibility[i][j][d] = 0;
      }
    }
  }
  return 0;
}

int write_up(void) {

 fprintf(outfile, "%.1f\t%.5f\t%.1f\t%.5f\n", correct_berf, berf_per_macro, correct_mberf, mrf_total_rep_per_macro);

 return 0;
}

int generate_demand() {

  for (int i=0; i<20000000; i++) {
    demand_list[i] = poisson3(demand_mean);
  }
  return 0;
}

int generate_prn() {

    for (int i=0; i<20000000; i++) {
      prn_list[i] = MRG32k3a3();
    }
    return 0;
  }

// int generate_demand(int sys_index, int constraint_index) {

//   for (int i=0; i<20000000; i++) {
//     demand_list[i] = poisson(demand_mean, sys_index, constraint_index);
//   }
//   return 0;
// }

// double generate_one_obs(int demand_index, int sys_index, int constraint_index) {

//   double total_cost =0;
//   double total_fail_prob = 0;
//   double num_stock_out_periods = 0;
//   double LittleS = system_value[sys_index][0];
//   int BigS = system_value[sys_index][1];
//   double current_level= BigS, next_level=0, Demand;
//   double Cost;

//   for(int i=0; i < 12; i++){
//     Cost = 0;
//     // Demand = poisson(demand_mean, sys_index, constraint_index);
//     Demand = demand_list[demand_index+i];

//     if( current_level < LittleS) {
//       next_level = BigS;
//       Cost = fixed_order_cost + order_cost * (BigS - current_level);
//     }
//     else next_level = current_level;

//     if( next_level - Demand >= 0) Cost += holding_cost * (next_level - Demand);
//     else  {
//       Cost += penalty_cost * (Demand - next_level);
//       total_fail_prob++;
//       num_stock_out_periods += 1;
//     }

//     current_level = next_level - Demand;
//     total_cost += Cost;
//   }

// //   printf("num stock outs: %.10f\n", num_stock_out_periods);

//   for (int j = 0; j < 1; j++) {
//     single_obs[0] = 0;
//     single_obs[1] = 0;
//     //printf("total cost: %.10f\n", total_cost);
//     if (total_cost > 1400) {
//       single_obs[0] = 1;
//     }

//     if (num_stock_out_periods >= 1) {
//       single_obs[1] = 1;
//     }
//   }

//   double prn = MRG32k3a_y(sys_index, constraint_index);
//   for (int j = 0; j < NumConstraint; j++) {
//     // printf("2nd obs: %.10f\n", single_obs[1]); inventory에서 dummies만들때 Numconstraint가 아니라 1로 들어가있는거 있음.. 두시간동안 찾았네.
//     for (int d = 0; d < NumThreshold; d++) {
//       //double rn = MRG32k3a();
//       dummies[j][d] = 0;
//       if (prn <= q[d][j]) {
//         dummies[j][d] = 1;
//       }
//     }
//   }
//   return 0;
// }

double generate_one_obs(int demand_index, int sys_index) {

  double total_cost =0;
  double total_fail_prob = 0;
  double num_stock_out_periods = 0;
  double LittleS = system_value[sys_index][0];
  int BigS = system_value[sys_index][1];
  double current_level= BigS, next_level=0, Demand;
  double Cost;

  for(int i=0; i < 12; i++){
    Cost = 0;

    if (CRN == 0){
        Demand = poisson2(demand_mean, sys_index);
    }
    
    else {
        Demand = demand_list[demand_index+i];
    }

    if( current_level < LittleS) {
      next_level = BigS;
      Cost = fixed_order_cost + order_cost * (BigS - current_level);
    }
    else next_level = current_level;

    if( next_level - Demand >= 0) Cost += holding_cost * (next_level - Demand);
    else  {
      Cost += penalty_cost * (Demand - next_level);
      total_fail_prob++;
      num_stock_out_periods += 1;
    }

    current_level = next_level - Demand;
    total_cost += Cost;
  }

//   printf("num stock outs: %.10f\n", num_stock_out_periods);

    single_obs[0] = 0;
    single_obs[1] = 0;
    //printf("total cost: %.10f\n", total_cost);
    if (total_cost > 1400) {
        single_obs[0] = 1;
    }

    if (num_stock_out_periods >= 1) {
        single_obs[1] = 1;
    }

  double prn = MRG32k3a_y2(sys_index);
  for (int j = 0; j < NumConstraint; j++) {
    for (int d = 0; d < NumThreshold; d++) {
      //double rn = MRG32k3a();
      dummies[j][d] = 0;
      if (prn <= q[d][j]) {
        dummies[j][d] = 1;
      }
    }
  }
  return 0;
}


double configuration(void) {

	// define system
    int s, S;
    int system_counter = 0;
    for (int j=0; j<11; j++) {
        s = 20 + 2 * j;
        for (int k=0; k<7; k++) {
            S = 40 + 10 * k;
            if (S >= s) {
                system_value[system_counter][0] = s;
                system_value[system_counter][1] = S;
                system_counter += 1;
            }
        }
    }

    if (NumConstraint == 1){
        q[0][0] = 0.01; 
        q[1][0] = 0.05; 
        q[2][0] = 0.1; 
        q[3][0] = 0.2;
    }

    else {
        q[0][0] = 0.01; 
        q[1][0] = 0.05; 
        q[2][0] = 0.1; 
        q[3][0] = 0.2;
        // q[4][0] = 0.3;
        // q[5][0] = 0.4;
        // q[6][0] = 0.5;

        q[0][1] = 0.01; 
        q[1][1] = 0.05; 
        q[2][1] = 0.1; 
        q[3][1] = 0.2;
        // q[4][1] = 0.3;
        // q[5][1] = 0.4;
        // q[6][1] = 0.5;
    }

	return 0;
}


double maxfn(double x, double y)
{
	if(x>y) return x;
	else return y;
}

double minfn(double x, double y)
{
	if(x<y) return x;
	else return y;
}

int main() {

    read_system_true_value();
    determine_true_feasibility();
    configuration();

    for (int i=0; i<NumSys; i++){
        printf("sys_true: %.10f\t%.10f\n", system_true_value[i][0], system_true_value[i][1]);
    }

    outfile = NULL;
    outfile = fopen("feasibiliy_MBeRF2_inventory","a");

    // for (int d = 0; d < NumThreshold; d++) {
    //     q[d][0] = 0.01 + 0.02 * (d);
    //     q[d][1] = 0.01 + 0.02 * (d);
    // }

    double alpha = 0.05;
    double beta = (1-pow(1-alpha, (double) 1/NumSys))/NumConstraint;    // independant systems
    double beta_prime = beta/2;

    for(int j = 0; j < NumConstraint; j++){
        H[j] = (log ((1/beta_prime) - 1))/(log (Theta));
        H[j] = std::ceil(H[j]);
        printf("H: %.1f\n", H[j]);
        for (int d=0; d<NumThreshold; d++) {
            printf("threshold %d %d : %.3f\n", d, j, q[d][j]);
        }
    }

    total_correct_berf = 0;
    total_correct_mberf = 0;
    total_match_decision = 0;

    double init_s10[NumSys];
    double init_s11[NumSys];
    double init_s12[NumSys];
    double init_s20[NumSys];
    double init_s21[NumSys];
    double init_s22[NumSys];

    double init_s10_y[NumSys];
    double init_s11_y[NumSys];
    double init_s12_y[NumSys];
    double init_s20_y[NumSys];
    double init_s21_y[NumSys];
    double init_s22_y[NumSys];

    berf_total = 0;
    mberf_total = 0;

    double matching_rep = 0;
    for (int l=0; l<NumMacro; l++) {
        
        std::random_device rd;
        std::mt19937 gen(rd());
        generate_demand();
        generate_prn();

        for (int i=0; i<NumSys; i++) {
            std::uniform_int_distribution<int> dist1(0, 4294967087);
            std::uniform_int_distribution<int> dist2(0, 4294944443);
            
            init_s10[i] = dist1(gen);
            init_s11[i] = dist1(gen);
            init_s12[i] = dist1(gen);
            init_s20[i] = dist2(gen);
            init_s21[i] = dist2(gen);
            init_s22[i] = dist2(gen);
        }

        // for (int i=0; i<NumSys; i++) {
        //         init_s10[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s11[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s12[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s20[i] = ((int) std::rand()) % ((4294944443 + 1));
        //         init_s21[i] = ((int) std::rand()) % ((4294944443 + 1));
        //         init_s22[i] = ((int) std::rand()) % ((4294944443 + 1));
        // }

        for (int i=0; i<NumSys; i++) {
                s10a[i] = init_s10[i];
                s11a[i] = init_s11[i];
                s12a[i] = init_s12[i];
                s20a[i] = init_s20[i];
                s21a[i] = init_s21[i];
                s22a[i] = init_s22[i];
        }

        for (int i=0; i<NumSys; i++) {
            std::uniform_int_distribution<int> dist1(0, 4294967087);
            std::uniform_int_distribution<int> dist2(0, 4294944443);
            
            init_s10_y[i] = dist1(gen);
            init_s11_y[i] = dist1(gen);
            init_s12_y[i] = dist1(gen);
            init_s20_y[i] = dist2(gen);
            init_s21_y[i] = dist2(gen);
            init_s22_y[i] = dist2(gen);
        }

        // for (int i=0; i<NumSys; i++) {
        //         init_s10_y[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s11_y[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s12_y[i] = ((int) std::rand()) % ((4294967087 + 1));
        //         init_s20_y[i] = ((int) std::rand()) % ((4294944443 + 1));
        //         init_s21_y[i] = ((int) std::rand()) % ((4294944443 + 1));
        //         init_s22_y[i] = ((int) std::rand()) % ((4294944443 + 1));
        // }

        for (int i=0; i<NumSys; i++) {
                s10_ya[i] = init_s10_y[i];
                s11_ya[i] = init_s11_y[i];
                s12_ya[i] = init_s12_y[i];
                s20_ya[i] = init_s20_y[i];
                s21_ya[i] = init_s21_y[i];
                s22_ya[i] = init_s22_y[i];
        }

//        printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10[0][0], s11[0][0], s12[0][0],s20[0][0], s21[0][0], s22[0][0]);

        for (int i=0; i<NumSys; i++) {

            int demand_index = 0;
            int prn_index = 0;

            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    BeRF_Z[i][j][d] = -2;
                    // printf("threshold %d %d : %.3f\n", d, j, q[d][j]);
                }
            }
        }

        // BeRF section
        correct_berf = 1;
        berf();
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    if (system_true_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * Theta))) {
                        if (BeRF_Z[i][j][d] == 1) correct_berf *= 1;
                        else correct_berf *= 0;
                    } else if (system_true_value[i][j] >= ((q[d][j] * Theta) / ((1 - q[d][j]) + q[d][j] * Theta))) {
                        if (BeRF_Z[i][j][d] == 0) correct_berf *= 1;
                        else correct_berf *= 0;
                    }
                }
            }
        }
        if (correct_berf == 1) total_correct_berf += 1;

        //printf("%d\t%d\t%d\t%d\n", RF_Z[0][0][0], RF_Z[0][0][1], RF_Z[0][0][2], RF_Z[0][0][3]);


        if(DOMBERF == 1){

            // MBeRF section
            correct_mberf = 1;

            // for (int p=0; p<NumPass; p++) {
            //     if (p = 0) {
            //         for (int j=0; j<NumConstraint; j++) {
            //             for (int d=0; d<NumThreshold; d++) {
            //                 if (d == 5 | d == 10 | d == 15 | d == 20 | d == 25) {
            //                     T_index[p][d][j] = 1;
            //                 }
            //                 else {
            //                     T_index[p][d][j] = 0;
            //                 }
            //             }
            //         }
            //     }
            //     if (p = 1) {
            //         for (int j=0; j<NumConstraint; j++) {
            //             for (int d=0; d<NumThreshold; d++) {
            //                 if (d == 0 | d == 1 | d == 2 | d == 3 | d == 4) {
            //                     T_index[p][d][j] = 1;
            //                 }
            //                 else {
            //                     T_index[p][d][j] = 0;
            //                 }
            //             }
            //         }
            //     }
            // }

            // first pass
            T_index[0][0][0] = 1; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 1; 
            T_index[0][0][1] = 1; T_index[0][1][1] = 0; T_index[0][2][1] = 0; T_index[0][3][1] = 1; 
            // T_index[0][4][0] = 1; T_index[0][5][0] = 1; T_index[0][6][0] = 1; T_index[0][7][0] = 1;

            // // second pass
            T_index[1][0][0] = 0; T_index[1][1][0] = 1; T_index[1][2][0] = 1; T_index[1][3][0] = 0;
            T_index[1][0][1] = 0; T_index[1][1][1] = 1; T_index[1][2][1] = 1; T_index[1][3][1] = 0;
            // T_index[1][4][0] = 0; T_index[1][5][0] = 0; T_index[1][6][0] = 0; T_index[1][7][0] = 0;

            // T_index[2][0][0] = 1; T_index[2][1][0] = 0; T_index[2][2][0] = 0; T_index[2][3][0] = 0;
            // T_index[2][0][1] = 1; T_index[2][1][1] = 0; T_index[2][2][1] = 0; T_index[2][3][1] = 0;

            // // first pass
            // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1; 
            // T_index[0][4][0] = 0; T_index[0][5][0] = 0; T_index[0][6][0] = 0; T_index[0][7][0] = 0;

            // // second pass
            // T_index[1][0][0] = 0; T_index[1][1][0] = 0; T_index[1][2][0] = 0; T_index[1][3][0] = 0;
            // T_index[1][4][0] = 1; T_index[1][5][0] = 1; T_index[1][6][0] = 1; T_index[1][7][0] = 1;

            // // first pass
            // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 0; T_index[0][3][0] = 0; 
            // T_index[0][4][0] = 0; T_index[0][5][0] = 0; T_index[0][6][0] = 1; T_index[0][7][0] = 1;

            // // second pass
            // T_index[1][0][0] = 0; T_index[1][1][0] = 0; T_index[1][2][0] = 1; T_index[1][3][0] = 1;
            // T_index[1][4][0] = 1; T_index[1][5][0] = 1; T_index[1][6][0] = 0; T_index[1][7][0] = 0;

            //all thresholds for 2 constraints. This is for 4 thresholds.
            // T_index[0][0][0] = 0; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 0; T_index[0][4][0] = 0; T_index[0][5][0] = 0; T_index[0][6][0] = 1;
            // T_index[0][0][1] = 0; T_index[0][1][1] = 0; T_index[0][2][1] = 0; T_index[0][3][1] = 0; T_index[0][4][1] = 0; T_index[0][5][1] = 0; T_index[0][6][1] = 1;

            // T_index[1][0][0] = 1; T_index[1][1][0] = 1; T_index[1][2][0] = 0; T_index[1][3][0] = 0; T_index[1][4][0] = 0; T_index[1][5][0] = 1; T_index[1][6][0] = 0;
            // T_index[1][0][1] = 1; T_index[1][1][1] = 1; T_index[1][2][1] = 0; T_index[1][3][1] = 0; T_index[1][4][1] = 0; T_index[1][5][1] = 1; T_index[1][6][1] = 0;

            // // T_index[2][0][0] = 0; T_index[2][1][0] = 0; T_index[2][2][0] = 0; T_index[2][3][0] = 0; T_index[2][4][0] = 1; T_index[2][5][0] = 0; T_index[2][6][0] = 0;
            // T_index[2][0][1] = 0; T_index[2][1][1] = 0; T_index[2][2][1] = 0; T_index[2][3][1] = 0; T_index[2][4][1] = 1; T_index[2][5][1] = 0; T_index[2][6][1] = 0;
            
            // T_index[3][0][0] = 0; T_index[3][1][0] = 0; T_index[3][2][0] = 0; T_index[3][3][0] = 1; T_index[3][4][0] = 0; T_index[3][5][0] = 0; T_index[3][6][0] = 0;
            // T_index[3][0][1] = 0; T_index[3][1][1] = 0; T_index[3][2][1] = 0; T_index[3][3][1] = 1; T_index[3][4][1] = 0; T_index[3][5][1] = 0; T_index[3][6][1] = 0;

            // T_index[4][0][0] = 0; T_index[4][1][0] = 0; T_index[4][2][0] = 1; T_index[4][3][0] = 0; T_index[4][4][0] = 0; T_index[4][5][0] = 0; T_index[4][6][0] = 0;
            // T_index[4][0][1] = 0; T_index[4][1][1] = 0; T_index[4][2][1] = 1; T_index[4][3][1] = 0; T_index[4][4][1] = 0; T_index[4][5][1] = 0; T_index[4][6][1] = 0;

            // T_index[5][0][0] = 0; T_index[5][1][0] = 1; T_index[5][2][0] = 0; T_index[5][3][0] = 0; T_index[5][4][0] = 0; T_index[5][5][0] = 0; T_index[5][6][0] = 0;
            // T_index[5][0][1] = 0; T_index[5][1][1] = 1; T_index[5][2][1] = 0; T_index[5][3][1] = 0; T_index[5][4][1] = 0; T_index[5][5][1] = 0; T_index[5][6][1] = 0;

            // T_index[6][0][0] = 1; T_index[6][1][0] = 0; T_index[6][2][0] = 0; T_index[6][3][0] = 0; T_index[6][4][0] = 0; T_index[6][5][0] = 0; T_index[6][6][0] = 0;
            // T_index[6][0][1] = 1; T_index[6][1][1] = 0; T_index[6][2][1] = 0; T_index[6][3][1] = 0; T_index[6][4][1] = 0; T_index[6][5][1] = 0; T_index[6][6][1] = 0;
            
            for (int i=0; i<NumSys; i++) {
                for (int j=0; j<NumConstraint; j++) {
                    for (int d=0; d<NumThreshold; d++) MBeRF_Z[i][j][d] = -2;
                }
            }

            for (int i=0; i<NumSys; i++) {
                    s10a[i] = init_s10[i];
                    s11a[i] = init_s11[i];
                    s12a[i] = init_s12[i];
                    s20a[i] = init_s20[i];
                    s21a[i] = init_s21[i];
                    s22a[i] = init_s22[i];
            }

            for (int i=0; i<NumSys; i++) {
                    s10_ya[i] = init_s10_y[i];
                    s11_ya[i] = init_s11_y[i];
                    s12_ya[i] = init_s12_y[i];
                    s20_ya[i] = init_s20_y[i];
                    s21_ya[i] = init_s21_y[i];
                    s22_ya[i] = init_s22_y[i];
            }
    //        printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10[0][0], s11[0][0], s12[0][0],s20[0][0], s21[0][0], s22[0][0]);

            mberf1(0);
        //    printf("%d\t%d\t%d\t%d\n", MRF_Z[0][0][0], MRF_Z[0][0][1], MRF_Z[0][0][2], MRF_Z[0][0][3]);
        //    printf("%d\t%d\t%d\t%d\n", MRF_Z[0][1][0], MRF_Z[0][1][1], MRF_Z[0][1][2], MRF_Z[0][1][3]);

            mberf2(1);
            // mberf2(2);
            // mberf2(3);
            // mberf2(4);
            // mberf2(5);
            // mberf2(6);
            
            // for (int i=0; i<NumSys; i++) {
            //     printf("system %d cons1 : %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, BeRF_Z[i][0][0], BeRF_Z[i][0][1], BeRF_Z[i][0][2], BeRF_Z[i][0][3]);
            //     printf("system %d cons1 : %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, MBeRF_Z[i][0][0], MBeRF_Z[i][0][1], MBeRF_Z[i][0][2], MBeRF_Z[i][0][3]);
            //     printf("system %d cons2 : %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, BeRF_Z[i][1][0], BeRF_Z[i][1][1], BeRF_Z[i][1][2], BeRF_Z[i][1][3]);
            //     printf("system %d cons2 : %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, MBeRF_Z[i][1][0], MBeRF_Z[i][1][1], MBeRF_Z[i][1][2], MBeRF_Z[i][1][3]);
            // }

            for (int i=0; i<NumSys; i++) {
                for (int j=0; j<NumConstraint; j++) {
                    for (int d=0; d<NumThreshold; d++) {
                    for (int p=0; p<NumPass; p++) {
                        if (T_index[p][d][j] == 1) {
                                if (system_true_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * Theta))) {
                                    if (MBeRF_Z[i][j][d] == 1) correct_mberf *= 1;
                                    else {correct_mberf *= 0;}
                                } else if (system_true_value[i][j] >= ((q[d][j] * Theta) / ((1 - q[d][j]) + q[d][j] * Theta))) {
                                    if (MBeRF_Z[i][j][d] == 0) correct_mberf *= 1;
                                    else correct_mberf *= 0;
                                }
                        }
                    }
                    }
                }
            }
            //printf("%.2f\n", correct_mrf);
            if (correct_mberf == 1) total_correct_mberf += 1;
        }

        double match_decision_indicator = 1;
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    for (int p=0; p<NumPass; p++) {
                        if (T_index[p][d][j] == 1) {
                            if (BeRF_Z[i][j][d] != MBeRF_Z[i][j][d]) {
                                match_decision_indicator *= 0;
                                // printf("system %d cons %d threshold %d %dth pass is not matched\n", i+1, j+1, d+1, p+1);
                            }
                        }
                    }
                }
            }
        }

        if (match_decision_indicator == 1) total_match_decision += 1;

        mrf_total_rep_per_macro = 0;
        for (int p=0; p<NumPass; p++) mrf_total_rep_per_macro += mberf_per_macro[p];
        if (mrf_total_rep_per_macro == berf_per_macro) matching_rep += 1;

        printf("%.1f\t%.5f\t%.1f\t%.5f\n", correct_berf, berf_per_macro, correct_mberf, mrf_total_rep_per_macro);
        
        write_up();
    
    }

    printf("BeRF: %.5f\n", total_correct_berf/NumMacro);
    for (int j=0; j<NumConstraint; j++) {
        printf("Constraint %d: %.5f\n", j+1, BeRF_total_obs[j]/NumMacro);
    }
    printf("BeRF total: %.5f\n", berf_total/NumMacro);

    if(DOMBERF==1){
        printf("MBeRF: %.5f\n", total_correct_mberf/NumMacro);
        for (int j=0; j<NumConstraint; j++) {
            printf("Constraint %d: %.5f\n", j+1, MBeRF_total_obs[j]/NumMacro);
        }

        for (int p=0; p<NumPass; p++) {
            printf("Pass %d: %.5f\n", p+1, MBeRF_rep_by_pass[p]/NumMacro);
        }

        printf("MBeRF total: %.5f\n", mberf_total/NumMacro);

        printf("Matching decision ratio: %.5f\n", total_match_decision/NumMacro);
        printf("Matching rep ratio: %.5f\n", matching_rep/NumMacro);
    }

    return 0;
}