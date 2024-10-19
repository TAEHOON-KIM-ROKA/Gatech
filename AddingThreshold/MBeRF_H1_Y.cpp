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

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define NumMacro 10000
#define NumSys	1
#define NumConstraint	2
#define NumThreshold	4
#define NumPass 2
#define Theta   1.5

// inputs for Generate R(0,1) by L'ecuyer (1997)
#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

double q[NumThreshold][NumConstraint];

double MRG32k3a(int sys_index, int constraint_index);  //Generate R(0,1) by L'ecuyer (1997)
// choices of seeds for Generate R(0,1) by L'ecuyer (1997)
double s10_y, s11_y, s12_y, s20_y, s21_y, s22_y;

double s10[NumSys][NumConstraint];
double s11[NumSys][NumConstraint];
double s12[NumSys][NumConstraint];
double s20[NumSys][NumConstraint];
double s21[NumSys][NumConstraint];
double s22[NumSys][NumConstraint];

double minfn(double x1, double x2);
double maxfn(double x, double y);
double Lower(double x, double y);
double Upper(double x, double y);

double normal(double rmean, double rvar, int sys_index, int constraint_index);
double configuration(void);
double generate_multiNormal(int numConstraint, int case_index);
int read_rand_seeds(void);
double generate_Bernoulli(int numConstraint, int sys_index, int constraint_index);

int berf(void);
int mberf1(int pass_index);
int mberf2(int pass_index);

double mean_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];

double prob[NumConstraint];
double H[NumConstraint];
double dummies[NumConstraint][NumThreshold];
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

double sumY[NumSys][NumConstraint];
double sumI[NumSys][NumConstraint][NumThreshold];
double num_obs[NumSys][NumConstraint];
double v_UB[NumSys][NumConstraint];
double v_LB[NumSys][NumConstraint];
int LAST[NumSys][NumConstraint];

double berf_total;
double mberf_total;
double berf_per_macro;
double mberf_per_macro[NumPass];
double mrf_total_rep_per_macro;

int lastthershold;

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

   		for (int j=0; j<NumConstraint; j++) {
   	        sumY[i][j] = 0;
   			for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] = 0;
            }
            surviveThreshold[j] = 0;
   		}

        for (int j = 0; j < NumConstraint; j++) {
            generate_Bernoulli(NumConstraint, i, j);
            sumY[i][j] += observations[i][j];
            for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] += dummies[j][d];
            }
            num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
            MBeRF_total_obs[j] += 1;
        }

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

                                if ((sumY[i][j]+H[j])/num_obs[i][j] <= (sumI[i][j][d]/num_obs[i][j])) {
       								MBeRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                else if ((sumY[i][j]-H[j])/num_obs[i][j] >= (sumI[i][j][d]/num_obs[i][j])) {
                                    MBeRF_Z[i][j][d] = 0;
       								ON_l[i][j][d] = 0;
       				                surviveThreshold[j] -= 1;
       							}
                                // printf("%.5f\t%.5f\t%.5f\t%d\n", v_LB[i][j], v_UB[i][j], (sumI[i][j][d]/num_obs[i][j]), d);
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
                // printf("MBeRF finished\n");
                break;
            }

   			for (int j=0; j<NumConstraint; j++) {

                // if (v_UB[i][j] > v_LB[i][j]) {
                if (ON[i][j] == 1) {
                    generate_Bernoulli(NumConstraint, i, j);
                    sumY[i][j] += observations[i][j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    MBeRF_total_obs[j] += 1;                
                // }
                }
            
   			}

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

                                        else if (LAST[i][j] == 1) {
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

            for (int j=0; j<NumConstraint; j++) {

                if (v_UB[i][j] > v_LB[i][j]) {
                    generate_Bernoulli(NumConstraint, i, j);
                    sumY[i][j] += observations[i][j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    MBeRF_total_obs[j] += 1;  
                }
   			}
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

                    // printf("%.5f\t%.5f\n", v_LB[i][j], v_UB[i][j]);

   					for (int d=0; d<NumThreshold; d++) {
    
                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                // printf("%.5f\t%.5f\t%.5f\n", v_LB[i][j], v_UB[i][j], (sumI[i][j][d]/num_obs[i][j]));

                                if (v_UB[i][j] <= q[d][j]) {
       								MBeRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                else if (v_LB[i][j] >= q[d][j]) {
                                    MBeRF_Z[i][j][d] = 0;
       								ON_l[i][j][d] = 0;
       				                surviveThreshold[j] -= 1;
       							}
                                // printf("%.5f\t%.5f\t%.5f\t%d\n", v_LB[i][j], v_UB[i][j], (sumI[i][j][d]/num_obs[i][j]), d);
                                if (surviveThreshold[j] == 1){
                                lastthershold = d;
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

        for (int j = 0; j < NumConstraint; j++) {
            generate_Bernoulli(NumConstraint, i, j);
            sumY[i][j] += observations[i][j];
            for (int d = 0; d < NumThreshold; d++) {
                sumI[i][j][d] += dummies[j][d];
            }
            num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
            BeRF_total_obs[j] += 1;
        }

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

                            if ((sumY[i][j]-H[j])/num_obs[i][j] >= (sumI[i][j][d]/num_obs[i][j])) {
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
                break;}

   			for (int j=0; j<NumConstraint; j++) {
                if (ON[i][j] == 1) {
                    generate_Bernoulli(NumConstraint, i, j);
                    sumY[i][j] += observations[i][j];
                    for (int d = 0; d < NumThreshold; d++) {
                        sumI[i][j][d] += dummies[j][d];
                    }
                    num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
                    BeRF_total_obs[j] += 1;
                }
   			}
            berf_total += 1;
            berf_per_macro += 1;
   		}

   	}

	return 0;
}

double MRG32k3a_y() //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11_y - a13n * s10_y;
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10_y = s11_y;   s11_y = s12_y;   s12_y = p1;

    // Component 2
    p2 = a21 * s22_y - a23n * s20_y;
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20_y = s21_y;   s21_y = s22_y;   s22_y = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
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

double generate_Bernoulli(int numConstraint, int sys_index, int constraint_index) {

    // if (NumConstraint < 1.5) {
        for (int i = 0; i < NumSys; i++) {
            double prn = MRG32k3a(sys_index, constraint_index);
            for (int j = 0; j < numConstraint; j++) {
                observations[i][j] = 0;
                if (MRG32k3a_y() <= prob[j]) {
                    observations[i][j] = 1;
                }
                //double prn = MRG32k3a();
                for (int d = 0; d < NumThreshold; d++) {
                    dummies[j][d] = 0;
                    if (prn <= q[d][j]) {
                        dummies[j][d] = 1;
                    }
                }
            }
            return 0;
        }
    // }

    // else if (NumConstraint > 1.5) {
    //     std::vector<std::vector<double>> C(numConstraint, std::vector<double>(numConstraint));

    //     for (int i = 0; i < numConstraint; i++) {
    //         for (int j = 0; j < numConstraint; j++) {
    //             C[i][j] = chol_matrix[0][i][j];
    //             //printf("C: % .2f\n", C[i][j]);
    //         }
    //     }

    //     for (int i = 0; i < NumSys; i++) {
    //         std::vector<double> std_normal(numConstraint);
    //         std::vector<double> correlated_normal(numConstraint);
    //         std::vector<int> successes(numConstraint);
    //         std::vector<double> x_value(numConstraint);

    //         for (int j = 0; j < numConstraint; j++){
    //             successes[j] = 0;
    //             boost::math::normal_distribution<> standard_normal;
    //             x_value[j] = boost::math::quantile(standard_normal, prob[j]);
    //         }

    //         for (int l = 0; l < numConstraint; l++) {
    //                 std_normal[l] = normal(0, 1, sys_index, constraint_index);
    //         }

    //         double prn = MRG32k3a(sys_index, constraint_index);
    //         for (int j = 0; j < numConstraint; j++) {

    //             correlated_normal[j] = 0;
    //             for (int m = 0; m < numConstraint; m++) {
    //                 correlated_normal[j] += C[j][m] * std_normal[m];
    //             }

    //             if (correlated_normal[j] < x_value[j]) {
    //                 successes[j]++;
    //             }

    //             //double prn = MRG32k3a();
    //             for (int d = 0; d < NumThreshold; d++) { 
    //                 dummies[j][d] = 0;
    //                 if (prn <= q[d][j]) {
    //                     dummies[j][d] = 1;
    //                 }
    //             }
    //         }

    //         for (int j = 0; j < numConstraint; j++) {
    //             observations[i][j] = static_cast<double>(successes[j]);
    //             //printf("obs: % .2f\n", observations[i][j]);
    //         }
    //         return 0;
    //     }
    // }
    return 0;
}

double configuration(void) {

    for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
			mean_value[i][j] = prob[j];
            // var_value[i][j] = 1;
		}
	}

    // Single system
    // for (int j=0; j<NumConstraint; j++) {
	// 	epsilon[j] = 1/sqrt(Nnot);
    //     //epsilon[j] = 0.1;
    // }

    // q[0][0] = Lower(prob, Theta*1.5); q[1][0] = Lower(prob, Theta); q[2][0] = Upper(prob, Theta); q[3][0] = Upper(prob, Theta*1.5);
    // q[0][1] = -3*epsilon[1]; q[1][1] = -epsilon[1]; q[2][1] = epsilon[1]; q[3][1] = 3*epsilon[1];

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

double Lower(double x, double y){
    // adjusting(x);
    // adjusting(y);
    double lower = (x)/(x + (1 - x)*y);
    // return adjusting(lower);
    return lower;
}

double Upper(double x, double y){
    // adjusting(x);
    // adjusting(y);
    double upper = (x*y)/(x*(y-1) + 1);
    // return adjusting(upper);
    return upper;
}

int write_up(void) {

 fprintf(outfile, "%.1f\t%.5f\t%.1f\t%.5f\n", berf_per_macro, mrf_total_rep_per_macro, lastthershold);

 return 0;
}

int main() {

    outfile = NULL;
    outfile = fopen("H1vsBeRF","a");

    prob[0] = 0.15;
    prob[1] = 0.15;

    q[0][0] = Lower(prob[0], Theta*1.5); q[1][0] = Lower(prob[0], Theta); q[2][0] = Upper(prob[0], Theta); q[3][0] = Upper(prob[0], Theta*1.5);
    q[0][1] = Lower(prob[1], Theta*1.5); q[1][1] = Lower(prob[1], Theta); q[2][1] = Upper(prob[1], Theta); q[3][1] = Upper(prob[1], Theta*1.5);

    double alpha = 0.05;
    double beta = (1-pow(1-alpha, (double) 1/NumSys))/NumConstraint;    // independant systems
    double beta_prime = beta/2;

    for(int j = 0; j < NumConstraint; j++){
        H[j] = (log ((1/beta_prime) - 1))/(log (Theta));
        H[j] = std::ceil(H[j]);
        printf("H: %.1f\n", H[j]);
        for(int d = 0; d < NumThreshold; d++){
            printf("threshold: %.8f\n", q[d][j]);
        }
    }

    double total_correct_berf = 0;
    double total_correct_mberf = 0;
    double total_match_decision = 0;

    double init_s10[NumSys][NumConstraint];
    double init_s11[NumSys][NumConstraint];
    double init_s12[NumSys][NumConstraint];
    double init_s20[NumSys][NumConstraint];
    double init_s21[NumSys][NumConstraint];
    double init_s22[NumSys][NumConstraint];

    double init_s10_y, init_s11_y, init_s12_y, init_s20_y, init_s21_y, init_s22_y;

    berf_total = 0;
    mberf_total = 0;

    double matching_rep = 0;
    for (int l=0; l<NumMacro; l++) {
        configuration();

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                init_s10[i][j] = ((int) std::rand()) % ((4294967087 + 1));
                init_s11[i][j] = ((int) std::rand()) % ((4294967087 + 1));
                init_s12[i][j] = ((int) std::rand()) % ((4294967087 + 1));
                init_s20[i][j] = ((int) std::rand()) % ((4294944443 + 1));
                init_s21[i][j] = ((int) std::rand()) % ((4294944443 + 1));
                init_s22[i][j] = ((int) std::rand()) % ((4294944443 + 1));
            }
        }

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                s10[i][j] = init_s10[i][j];
                s11[i][j] = init_s11[i][j];
                s12[i][j] = init_s12[i][j];
                s20[i][j] = init_s20[i][j];
                s21[i][j] = init_s21[i][j];
                s22[i][j] = init_s22[i][j];
            }
        }


        init_s10_y = ((int) std::rand()) % ((4294967087 + 1));
        init_s11_y = ((int) std::rand()) % ((4294967087 + 1));
        init_s12_y = ((int) std::rand()) % ((4294967087 + 1));
        init_s20_y = ((int) std::rand()) % ((4294967087 + 1));
        init_s21_y = ((int) std::rand()) % ((4294967087 + 1));
        init_s22_y = ((int) std::rand()) % ((4294967087 + 1));

        s10_y = init_s10_y, s11_y = init_s11_y, s12_y = init_s12_y, s20_y = init_s20_y, s21_y = init_s21_y, s22_y = init_s22_y;

        // if (l == 0) printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10_y, s11_y, s12_y, s20_y, s21_y, s22_y);

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) BeRF_Z[i][j][d] = -2;
            }
        }

        // BeRF section
        double correct_berf = 1;
        berf();
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    if (mean_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * Theta))) {
                        if (BeRF_Z[i][j][d] == 1) correct_berf *= 1;
                        else correct_berf *= 0;
                    } else if (mean_value[i][j] >= ((q[d][j] * Theta) / ((1 - q[d][j]) + q[d][j] * Theta))) {
                        if (BeRF_Z[i][j][d] == 0) correct_berf *= 1;
                        else correct_berf *= 0;
                    }
                }
            }
        }
        if (correct_berf == 1) total_correct_berf += 1;

        //printf("%d\t%d\t%d\t%d\n", RF_Z[0][0][0], RF_Z[0][0][1], RF_Z[0][0][2], RF_Z[0][0][3]);

        // MBeRF section
        double correct_mberf = 1;

        // // first pass
        // T_index[0][0][0] = 1; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 1;  // T_index[NumPass][NumThreshold][NumConstraint]
        // T_index[0][0][1] = 1; T_index[0][1][1] = 0; T_index[0][2][1] = 0; T_index[0][3][1] = 1;
        // //T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1;

        // // second pass
        // T_index[1][0][0] = 0; T_index[1][1][0] = 1; T_index[1][2][0] = 1; T_index[1][3][0] = 0;
        // T_index[1][0][1] = 0; T_index[1][1][1] = 1; T_index[1][2][1] = 1; T_index[1][3][1] = 0;
        // //T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        // first pass
        T_index[0][0][0] = 0; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 0;  // T_index[NumPass][NumThreshold][NumConstraint]
        T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        //T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1;

        // second pass
        T_index[1][0][0] = 1; T_index[1][1][0] = 0; T_index[1][2][0] = 0; T_index[1][3][0] = 1;
        T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;
        //T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        // // first pass
        // T_index[0][0][0] = 0; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 0;  // T_index[NumPass][NumThreshold][NumConstraint]
        // T_index[0][0][1] = 1; T_index[0][1][1] = 0; T_index[0][2][1] = 0; T_index[0][3][1] = 1;
        // //T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1;

        // // second pass
        // T_index[1][0][0] = 1; T_index[1][1][0] = 0; T_index[1][2][0] = 0; T_index[1][3][0] = 1;
        // T_index[1][0][1] = 0; T_index[1][1][1] = 1; T_index[1][2][1] = 1; T_index[1][3][1] = 0;
        // //T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        //         // first pass
        // T_index[0][0][0] = 1; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 1;  // T_index[NumPass][NumThreshold][NumConstraint]
        // T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // //T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1;

        // // second pass
        // T_index[1][0][0] = 0; T_index[1][1][0] = 1; T_index[1][2][0] = 1; T_index[1][3][0] = 0;
        // T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;
        // // T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        // // first pass
        // T_index[0][0][0] = 1; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 1;  // T_index[NumPass][NumThreshold][NumConstraint]
        // T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // //T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;
        // // T_index[0][0][0] = 1; T_index[0][1][0] = 1; T_index[0][2][0] = 1; T_index[0][3][0] = 1;

        // // second pass
        // T_index[1][0][0] = 0; T_index[1][1][0] = 1; T_index[1][2][0] = 1; T_index[1][3][0] = 0;
        // T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;
        // //T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) MBeRF_Z[i][j][d] = -2;
            }
        }

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                s10[i][j] = init_s10[i][j];
                s11[i][j] = init_s11[i][j];
                s12[i][j] = init_s12[i][j];
                s20[i][j] = init_s20[i][j];
                s21[i][j] = init_s21[i][j];
                s22[i][j] = init_s22[i][j];
            }
        }
        s10_y = init_s10_y, s11_y = init_s11_y, s12_y = init_s12_y, s20_y = init_s20_y, s21_y = init_s21_y, s22_y = init_s22_y;

        // if (l == 0) printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10[0][0], s11[0][0], s12[0][0],s20[0][0], s21[0][0], s22[0][0]);

        mberf1(0);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][0][0], MRF_Z[0][0][1], MRF_Z[0][0][2], MRF_Z[0][0][3]);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][1][0], MRF_Z[0][1][1], MRF_Z[0][1][2], MRF_Z[0][1][3]);

        mberf2(1);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][0][0], MRF_Z[0][0][1], MRF_Z[0][0][2], MRF_Z[0][0][3]);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][1][0], MRF_Z[0][1][1], MRF_Z[0][1][2], MRF_Z[0][1][3]);

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                //    for (int p=0; p<NumPass; p++) {
                //        if (T_index[p][d][j] == 1) {
                            if (mean_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * Theta))) {
                                if (MBeRF_Z[i][j][d] == 1) correct_mberf *= 1;
                                else {correct_mberf *= 0;}
                            } else if (mean_value[i][j] >= ((q[d][j] * Theta) / ((1 - q[d][j]) + q[d][j] * Theta))) {
                                if (MBeRF_Z[i][j][d] == 0) correct_mberf *= 1;
                                else correct_mberf *= 0;
                            }
                //        }
                //    }
                }
            }
        }
        //printf("%.2f\n", correct_mrf);
        if (correct_mberf == 1) total_correct_mberf += 1;

        double match_decision_indicator = 1;
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    for (int p=0; p<NumPass; p++) {
                        if (T_index[p][d][j] == 1) {
                            if (BeRF_Z[i][j][d] != MBeRF_Z[i][j][d]) match_decision_indicator *= 0;
                        }
                    }
                }
            }
        }

        if (match_decision_indicator == 1) total_match_decision += 1;

        mrf_total_rep_per_macro = 0;
        for (int p=0; p<NumPass; p++) mrf_total_rep_per_macro += mberf_per_macro[p];
        // printf("%.5f\t%.5f\t%d\n", berf_per_macro, mrf_total_rep_per_macro, lastthershold);
        write_up();
        if (mrf_total_rep_per_macro == berf_per_macro) matching_rep += 1;

    }

    printf("BeRF: %.5f\n", total_correct_berf/NumMacro);
    for (int j=0; j<NumConstraint; j++) {
        printf("Constraint %d: %.5f\n", j+1, BeRF_total_obs[j]/NumMacro);
    }
    printf("BeRF total: %.5f\n", berf_total/NumMacro);

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

    printf("Estimated p probability: %.5f\t%.5f\t%.5f\t%.5f\n", sumI[0][0][0]/num_obs[0][0], sumI[0][0][1]/num_obs[0][0],
                                                                sumI[0][0][2]/num_obs[0][0], sumI[0][0][3]/num_obs[0][0]);

    printf("Estimated num obs: %.5f\n", num_obs[0][0]);

    return 0;
}