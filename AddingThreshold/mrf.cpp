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

using namespace std;

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define Nnot	20
#define NumMacro 10000
#define NumSys	1
#define NumConstraint	2
#define NumThreshold	4
#define NumPass 2

// inputs for Generate R(0,1) by L'ecuyer (1997)
#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0


double MRG32k3a(int sys_index, int constraint_index);  //Generate R(0,1) by L'ecuyer (1997)
// choices of seeds for Generate R(0,1) by L'ecuyer (1997)
//double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;

double s10[NumSys][NumConstraint];
double s11[NumSys][NumConstraint];
double s12[NumSys][NumConstraint];
double s20[NumSys][NumConstraint];
double s21[NumSys][NumConstraint];
double s22[NumSys][NumConstraint];

double minfn(double x1, double x2);
double maxfn(double x, double y);

double normal(double rmean, double rvar, int sys_index, int constraint_index);
double configuration(void);
double generate_multiNormal(int numConstraint, int case_index);
int read_rand_seeds(void);

int rf(void);
int mrf(int pass_index);
int mrf2(int pass_index);

double mean_value[NumSys][NumConstraint];
double var_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];
double q[NumThreshold][NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int RF_Z[NumSys][NumConstraint][NumThreshold];
int MRF_Z[NumSys][NumConstraint][NumThreshold];
double RF_total_obs[NumConstraint];
double MRF_total_obs[NumConstraint];
double MRF_rep_by_pass[NumPass];
int T_index[NumPass][NumThreshold][NumConstraint];

double sumY[NumSys][NumConstraint];
double sum_squareY[NumSys][NumConstraint];
double num_obs[NumSys][NumConstraint];
double R[NumSys][NumConstraint];
double Sil2[NumSys][NumConstraint];
double v_UB[NumSys][NumConstraint];
double v_LB[NumSys][NumConstraint];
int LAST[NumSys][NumConstraint];

double eta[NumConstraint];

double rf_total;
double mrf_total;
double rf_per_macro;
double mrf_per_macro[NumPass];

int main() {

    double alpha = 0.05;
    double beta = (1-pow(1-alpha, (double) 1/NumSys))/NumConstraint;    // independant systems
    double beta_prime = beta/2;

    for (int j=0; j<NumConstraint; j++) eta[j] = 0.5*(pow(2*beta_prime, -2/((double) Nnot-1))-1);
    printf("%.5f\n", eta[0]);

    double total_correct_rf = 0;
    double total_correct_mrf = 0;
    double total_match_decision = 0;

    double init_s10[NumSys][NumConstraint];
    double init_s11[NumSys][NumConstraint];
    double init_s12[NumSys][NumConstraint];
    double init_s20[NumSys][NumConstraint];
    double init_s21[NumSys][NumConstraint];
    double init_s22[NumSys][NumConstraint];

    rf_total = 0;
    mrf_total = 0;

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

//        printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10[0][0], s11[0][0], s12[0][0],s20[0][0], s21[0][0], s22[0][0]);

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) RF_Z[i][j][d] = -2;
            }
        }

        // RF section
        double correct_rf = 1;
        rf();
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    if (mean_value[i][j] <= q[d][j]-epsilon[j]) {
                        if (RF_Z[i][j][d] == 1) correct_rf *= 1;
                        else correct_rf *= 0;
                    } else if (mean_value[i][j] >= q[d][j]+epsilon[j]) {
                        if (RF_Z[i][j][d] == 0) correct_rf *= 1;
                        else correct_rf *= 0;
                    }
                }
            }
        }
        if (correct_rf == 1) total_correct_rf += 1;

        //printf("%d\t%d\t%d\t%d\n", RF_Z[0][0][0], RF_Z[0][0][1], RF_Z[0][0][2], RF_Z[0][0][3]);

        // MRF section
        double correct_mrf = 1;

        // first pass
        T_index[0][0][0] = 1; T_index[0][1][0] = 0; T_index[0][2][0] = 0; T_index[0][3][0] = 1;
        T_index[0][0][1] = 0; T_index[0][1][1] = 1; T_index[0][2][1] = 1; T_index[0][3][1] = 0;

        // second pass
        T_index[1][0][0] = 0; T_index[1][1][0] = 1; T_index[1][2][0] = 1; T_index[1][3][0] = 0;
        T_index[1][0][1] = 1; T_index[1][1][1] = 0; T_index[1][2][1] = 0; T_index[1][3][1] = 1;

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) MRF_Z[i][j][d] = -2;
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
//        printf("%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", s10[0][0], s11[0][0], s12[0][0],s20[0][0], s21[0][0], s22[0][0]);

        mrf(0);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][0][0], MRF_Z[0][0][1], MRF_Z[0][0][2], MRF_Z[0][0][3]);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][1][0], MRF_Z[0][1][1], MRF_Z[0][1][2], MRF_Z[0][1][3]);

        mrf2(1);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][0][0], MRF_Z[0][0][1], MRF_Z[0][0][2], MRF_Z[0][0][3]);
//        printf("%d\t%d\t%d\t%d\n", MRF_Z[0][1][0], MRF_Z[0][1][1], MRF_Z[0][1][2], MRF_Z[0][1][3]);

        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                //    for (int p=0; p<NumPass; p++) {
                //        if (T_index[p][d][j] == 1) {
                            if (mean_value[i][j] <= q[d][j]-epsilon[j]) {
                                if (MRF_Z[i][j][d] == 1) correct_mrf *= 1;
                                else {correct_mrf *= 0;}
                            } else if (mean_value[i][j] >= q[d][j]+epsilon[j]) {
                                if (MRF_Z[i][j][d] == 0) correct_mrf *= 1;
                                else correct_mrf *= 0;
                            }
                //        }
                //    }
                }
            }
        }
        //printf("%.2f\n", correct_mrf);
        if (correct_mrf == 1) total_correct_mrf += 1;

        double match_decision_indicator = 1;
        for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                    for (int p=0; p<NumPass; p++) {
                        if (T_index[p][d][j] == 1) {
                            if (RF_Z[i][j][d] != MRF_Z[i][j][d]) match_decision_indicator *= 0;
                        }
                    }
                }
            }
        }

        if (match_decision_indicator == 1) total_match_decision += 1;

        double mrf_total_rep_per_macro = 0;
        for (int p=0; p<NumPass; p++) mrf_total_rep_per_macro += mrf_per_macro[p];
        if (mrf_total_rep_per_macro == rf_per_macro) matching_rep += 1;

    }

    printf("RF: %.5f\n", total_correct_rf/NumMacro);
    for (int j=0; j<NumConstraint; j++) {
        printf("Constraint %d: %.5f\n", j+1, RF_total_obs[j]/NumMacro);
    }
    printf("RF total: %.5f\n", rf_total/NumMacro);

    printf("MRF: %.5f\n", total_correct_mrf/NumMacro);
    for (int j=0; j<NumConstraint; j++) {
        printf("Constraint %d: %.5f\n", j+1, MRF_total_obs[j]/NumMacro);
    }

    for (int p=0; p<NumPass; p++) {
        printf("Pass %d: %.5f\n", p+1, MRF_rep_by_pass[p]/NumMacro);
    }

    printf("MRF total: %.5f\n", mrf_total/NumMacro);

    printf("Matching decision ratio: %.5f\n", total_match_decision/NumMacro);
    printf("Matching rep ratio: %.5f\n", matching_rep/NumMacro);

    return 0;
}

int mrf(int pass_index) {

    mrf_per_macro[pass_index] = 0;

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
            Sil2[i][j] = 0;
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
   			sum_squareY[i][j] = 0;
            surviveThreshold[j] = 0;
   		}

   		for (int n=0; n<Nnot; n++) {

   			for (int j=0; j<NumConstraint; j++) {
                observations[i][j] = normal(mean_value[i][j], var_value[i][j], i, j);
   				sumY[i][j] += observations[i][j];
   				sum_squareY[i][j] += observations[i][j] * observations[i][j];
                num_obs[i][j] += 1;
                MRF_total_obs[j] += 1;
   			}

            mrf_total += 1;
            MRF_rep_by_pass[pass_index] += 1;
            mrf_per_macro[pass_index] += 1;

   		}

        // for same number of thresholds per constriant
        for (int j=0; j<NumConstraint; j++) {
            for (int d=0; d<NumThreshold; d++) surviveThreshold[j] += T_index[pass_index][d][j];
        }

   		// find continuation region
   		for (int j=0; j<NumConstraint; j++) {
   			Sil2[i][j] = (sum_squareY[i][j]/(Nnot-1)) - (sumY[i][j]/(Nnot-1))*(sumY[i][j]/Nnot);
   			R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
   		}

   		while (surviveConstraint != 0) {

   			for (int j=0; j<NumConstraint; j++) {

   				if ((v_UB[i][j] > v_LB[i][j])) {

                    if ((sumY[i][j]+R[i][j])/num_obs[i][j] < v_UB[i][j]) LAST[i][j] = 0;   // UB is 0
                    if ((sumY[i][j]-R[i][j])/num_obs[i][j] > v_LB[i][j]) LAST[i][j] = 1;   // LB is 1

                    v_UB[i][j] = minfn(v_UB[i][j], (sumY[i][j]+R[i][j])/num_obs[i][j]);
                    v_LB[i][j] = maxfn(v_LB[i][j], (sumY[i][j]-R[i][j])/num_obs[i][j]);

   					for (int d=0; d<NumThreshold; d++) {

                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                if (v_UB[i][j] <= q[d][j]) {
       								MRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                if (v_LB[i][j] >= q[d][j]) {
                                    MRF_Z[i][j][d] = 0;
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
   				}
   			}

   			if (surviveConstraint == 0) break;

   			for (int j=0; j<NumConstraint; j++) {

                if (v_UB[i][j] > v_LB[i][j]) {
                    observations[i][j] = normal(mean_value[i][j], var_value[i][j], i, j);
       				sumY[i][j] += observations[i][j];
                    num_obs[i][j] += 1;
                    MRF_total_obs[j] += 1;
                    R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
                }
   			}

            mrf_total += 1;
            MRF_rep_by_pass[pass_index] += 1;
            mrf_per_macro[pass_index] += 1;
   		}
   	}

    return 0;
}

int mrf2(int pass_index) {

    mrf_per_macro[pass_index] = 0;

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
                if (T_index[i][d][j] == 1) {
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

                                if (((v_UB[i][j] <= q[d][j]) & (q[d][j] < v_LB[i][j])) |
                                    ((v_UB[i][j] < q[d][j]) & (q[d][j] <= v_LB[i][j]))) {

                                        if (LAST[i][j] == 0) {
                                            MRF_Z[i][j][d] = 0;
                                            ON_l[i][j][d] = 0;
                                            surviveThreshold[j] -= 1;
                                        }

                                        if (LAST[i][j] == 1) {
                                            MRF_Z[i][j][d] = 1;
                                            ON_l[i][j][d] = 0;
                                            surviveThreshold[j] -= 1;
                                        }
                                } else if ((v_UB[i][j] <= q[d][j]) & (v_LB[i][j] <= q[d][j])) {
                                    MRF_Z[i][j][d] = 1;
                                    ON_l[i][j][d] = 0;
                                    surviveThreshold[j] -= 1;
                                } else if ((v_LB[i][j] >= q[d][j]) & (v_UB[i][j] >= q[d][j])) {
                            //        if (j == 0) printf("%.5f\t%.5f\t%d\n", v_LB[i][j], v_UB[i][j], d);
                                    MRF_Z[i][j][d] = 0;
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
            if (surviveConstraint == 0) break;

            for (int j=0; j<NumConstraint; j++) {

                if (v_UB[i][j] > v_LB[i][j]) {
                    observations[i][j] = normal(mean_value[i][j], var_value[i][j], i, j);
       				sumY[i][j] += observations[i][j];
                    num_obs[i][j] += 1;
                    MRF_total_obs[j] += 1;
                    R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
                }
   			}
            mrf_total += 1;
            MRF_rep_by_pass[pass_index] += 1;
            mrf_per_macro[pass_index] += 1;

            for (int j=0; j<NumConstraint; j++) {

    //            if (j == 0) printf("%.5f\t%.5f\n", v_LB[i][j], v_UB[i][j]);

   				if ((ON[i][j] == 1) & (v_UB[i][j] > v_LB[i][j])) {

                    if ((sumY[i][j]+R[i][j])/num_obs[i][j] < v_UB[i][j]) LAST[i][j] = 0;   // UB is 0
                    if ((sumY[i][j]-R[i][j])/num_obs[i][j] > v_LB[i][j]) LAST[i][j] = 1;   // LB is 1

                    v_UB[i][j] = minfn(v_UB[i][j], (sumY[i][j]+R[i][j])/num_obs[i][j]);
                    v_LB[i][j] = maxfn(v_LB[i][j], (sumY[i][j]-R[i][j])/num_obs[i][j]);

   					for (int d=0; d<NumThreshold; d++) {

                        if (T_index[pass_index][d][j] == 1) {

                            if (ON_l[i][j][d] == 1) {

                                if (v_UB[i][j] <= q[d][j]) {
       								MRF_Z[i][j][d] = 1;
       								ON_l[i][j][d] = 0;
       								surviveThreshold[j] -= 1;
                                }

                                if (v_LB[i][j] >= q[d][j]) {
                                    MRF_Z[i][j][d] = 0;
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

   			if (surviveConstraint == 0) break;
        }
    }

    return 0;
}

int rf(void) {

    rf_per_macro = 0;

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
            Sil2[i][j] = 0;
            num_obs[i][j] = 0;
   		}
   	}

   	for (int i=0; i<NumSys; i++) {

   	    // generate initial samples
        int surviveConstraint = NumConstraint;
        int surviveThreshold[NumConstraint];

   		for (int j=0; j<NumConstraint; j++) {
   	        sumY[i][j] = 0;
   			sum_squareY[i][j] = 0;
   		}

   		for (int n=0; n<Nnot; n++) {

   			for (int j=0; j<NumConstraint; j++) {
                observations[i][j] = normal(mean_value[i][j], var_value[i][j], i, j);
   				sumY[i][j] += observations[i][j];
   				sum_squareY[i][j] += observations[i][j] * observations[i][j];
                num_obs[i][j] += 1;
                RF_total_obs[j] += 1;
   			}
            rf_total += 1;
            rf_per_macro += 1;

   		}

        // for same number of thresholds per constriant
        for (int j=0; j<NumConstraint; j++) surviveThreshold[j] = NumThreshold;

   		// find continuation region
   		for (int j=0; j<NumConstraint; j++) {
   			Sil2[i][j] = (sum_squareY[i][j]/(Nnot-1)) - (sumY[i][j]/(Nnot-1))*(sumY[i][j]/Nnot);
   			R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
   		}

   		while (surviveConstraint != 0) {

   			for (int j=0; j<NumConstraint; j++) {

   				if (ON[i][j] == 1) {

   					for (int d=0; d<NumThreshold; d++) {

   						if (ON_l[i][j][d] == 1) {

                            if ((sumY[i][j]+R[i][j])/num_obs[i][j] <= q[d][j]) {
   								RF_Z[i][j][d] = 1;
   								ON_l[i][j][d] = 0;
   								surviveThreshold[j] -= 1;
                            }

                            if ((sumY[i][j]-R[i][j])/num_obs[i][j] >= q[d][j]) {
                                RF_Z[i][j][d] = 0;
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

   			for (int j=0; j<NumConstraint; j++) {
                if (ON[i][j] == 1) {
                    observations[i][j] = normal(mean_value[i][j], var_value[i][j], i, j);
       				sumY[i][j] += observations[i][j];
                    num_obs[i][j] += 1;
                    RF_total_obs[j] += 1;
                    R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
                }
   			}
            rf_total += 1;
            rf_per_macro += 1;
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

double configuration(void) {

    for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
			mean_value[i][j] = 0;
            var_value[i][j] = 1;
		}
	}

    // Single system
    for (int j=0; j<NumConstraint; j++) {
		epsilon[j] = 1/sqrt(Nnot);
        //epsilon[j] = 0.1;
    }

    q[0][0] = -3*epsilon[0]; q[1][0] = -epsilon[0]; q[2][0] = epsilon[0]; q[3][0] = 3*epsilon[0];
    q[0][1] = -3*epsilon[1]; q[1][1] = -epsilon[1]; q[2][1] = epsilon[1]; q[3][1] = 3*epsilon[1];

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
