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
#define alpha 0.05
#define Nnot	20
#define NumMacro 100
#define NumSys	77
#define NumConstraint	2
#define NumThreshold	4
#define Num_s 20
#define Num_S 20
#define NumBatch    32
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
//double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;
//double  s10 = 43, s11 =54, s12 =65, s20 =43, s21 =54, s22 =65;
//double  s10 = 4321111, s11 =1115432, s12 =1116543, s20 =4321111, s21 =1115432, s22 =6543111;
//double  s10 = 43221, s11 =54332, s12 =65443, s20 =43321, s21 =54532, s22 =61543;
//double  s10 = 1010, s11 =10, s12 =101, s20 =2001, s21 = 202, s22 = 202;
//double  s10 = 4321, s11 =5432, s12 =6543, s20 =4321, s21 =5432, s22 =6543;
//double  s10 = 11, s11 =12, s12 =13, s20 =21, s21 =22, s22 =23;
double  s10 = 20202020, s11 =10101010, s12 =30303030, s20 =50505050, s21 =70707070, s22 =90909090;
//double  s10 = 151234, s11 =245245, s12 =1234134, s20 =34563, s21 =12341, s22 =234524;
//double  s10 = 23411, s11 = 24231, s12 = 41231, s20 = 31312, s21 = 42312, s22 = 32312;

double minfn(double x1, double x2);
double maxfn(double x, double y);

double normal(double rmean, double rvar);
double poisson(double pmean);
double configuration(void);
double generate_one_obs(int system_index, int demand_index);
int generate_demand(void);
int write_up(void);
int read_system_true_value(void);
int determine_true_feasibility(void);

double q[NumConstraint][NumThreshold];
double qL[NumConstraint][NumThreshold];
double qU[NumConstraint][NumThreshold];
double newq[NumConstraint][NumThreshold];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];

int system_value[NumSys][2];
double system_true_value[NumSys][2];
int true_feasibility[NumSys][NumConstraint][NumThreshold];
double single_obs[NumConstraint];
double single_sub_obs[NumConstraint];
double demand_list[2000000];

double demand_mean = 25;
double order_cost = 3;
double fixed_order_cost = 32;
double holding_cost = 1;
double penalty_cost = 5;

double total_obs;
double final_cd;
double batch_size = NumBatch;

FILE *outfile;

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
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20 = s21;   s21 = s22;   s22 = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double normal(double rmean, double rvar)
/* return normal random variable with mean rmean and variance rvar. */
// this is modified for Fully Sequential Procedure with CRN
{
	double V1 = 0, V2 = 0, W = 2, Y = 0, X1 = 0;
	do {
		V1 = 2 * MRG32k3a() - 1;
		V2 = 2 * MRG32k3a() - 1;
     	W = pow(V1,2) + pow(V2,2);

	} while (W > 1);

	Y = sqrt( (-2.00 * log(W))/W );
	X1 = rmean + sqrt(rvar) * V1 * Y;
	return X1;
}

double poisson(double lam) {
    double a, b;
    int i;
    a=exp(-lam);
    b=1;
    i=0;

    while(1) {
        b=b*MRG32k3a();
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

  return 0;
}

int determine_true_feasibility(void) {

  for (int i=0; i<NumSys; i++) {
    for (int j=0; j<NumConstraint; j++) {
      for (int d=0; d<NumThreshold; d++) {
        if (system_true_value[i][j] <= q[j][d]) true_feasibility[i][j][d] = 1;
        else true_feasibility[i][j][d] = 0;
      }
    }
  }
  return 0;
}

int write_up(void) {

 fprintf(outfile, "%.1f\t%.5f\n", total_obs, final_cd);

 return 0;
}

int generate_demand() {

  for (int i=0; i<2000000; i++) {
    demand_list[i] = poisson(demand_mean);
  }
  return 0;
}

double generate_one_obs(int system_index, int demand_index) {

  	double total_fail_prob = 0;
    double num_stock_out_periods = 0;
  	double LittleS = system_value[system_index][0];
  	int BigS = system_value[system_index][1];
  	double current_level= BigS, next_level=0, Demand;
  	double Cost;
  	double total_ber = 0;
    double total_stockout = 0;
  	//for (int b=0; b<batch_size; b++) {
    double total_cost =0;
    for(int t=0; t < 12; t++) {
      Cost = 0;
      Demand = poisson(demand_mean);
      //Demand = demand_list[demand_index+j];

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

    	//if (system_index == 18) printf("%.4f\t%.4f\t%f\n", total_cost, Cost, total_ber);

      single_sub_obs[0] = 0;
      single_sub_obs[1] = 0;

    	if (total_cost > 1400) {
        single_sub_obs[0] = 1;
      }

      if (num_stock_out_periods >= 1) {
        single_sub_obs[1] = 1;
      }
  	//}

  	//single_obs[0] = total_ber/batch_size;
  	//if (system_index == 18) printf("%.2f\t%.4f\n", yearly_cost_greater_than1400, batch_size);
  	//if (system_index == 18) 
  	//	printf("%.4f\t%f\t%d\t%f\t%f\n", single_obs[0], LittleS, BigS, total_ber, batch_size);


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

    for (int i=0; i<NumSys; i++) {
        for (int j=0; j<NumConstraint; j++) {
            ON[i][j] = 1;

            for (int d=0; d<NumThreshold; d++) {
                ON_l[i][j][d] = 1;
            }
        }
    }

    
    if (Theta == 1.5){
      for (int j=0; j<NumConstraint; j++) {  
        epsilon[j] = 0.004118;  // odd ratio = 1.5
      }
    }
    else if (Theta == 1.2){
      for (int j=0; j<NumConstraint; j++) {
        epsilon[j] = 0.001814;  // odd ratio = 1.2
      }
    }
  // set threshold value. The smallest first.
    q[0][0] = 0.01;
    q[0][1] = 0.05;
    q[0][2] = 0.1;
    q[0][3] = 0.2;

    q[1][0] = 0.01;
    q[1][1] = 0.05;
    q[1][2] = 0.1;
    q[1][3] = 0.2;

    for (int j = 0; j < NumConstraint; j++) {
        double UB = q[j][0]*Theta/(q[j][0]*Theta + (1-q[j][0]));
        double LB = q[j][0]/(q[j][0] + (1-q[j][0])*Theta);
        epsilon[j] = (UB-LB)/2;
        //printf("epsilon: %.4f\n", epsilon[j]);
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

int main()
{
    read_system_true_value();
    determine_true_feasibility();

    for (int i=0; i<NumSys; i++){
        printf("sys %d_true: %.10f\t%.10f\n", i, system_true_value[i][0], system_true_value[i][1]);
    }

    outfile = NULL;
    outfile = fopen("RF_2const_32_1.5_mock.out","a");

    double eta[NumConstraint];
    double beta = (1-pow(1-alpha, (double) 1/NumSys));
    //double beta = (double) alpha/NumSys;
    //printf("%.4f\n", beta);
    // if (NumThreshold == 1) beta = beta;
    beta = beta/(2*NumConstraint);
    for (int j=0; j<NumConstraint; j++) eta[j] = 0.5*( pow( 2*beta,(double) -2/(Nnot-1)) - 1);

    printf("%.4f\n", eta[0]);
    printf("%.4f\n", eta[1]);

	  double total_single_obs[NumConstraint];

    for (int l=0; l<NumMacro; l++) {

        configuration();
        total_obs = 0;
        final_cd = 1;

        double num_obs[NumSys][NumConstraint];
        double R[NumSys][NumConstraint];
        double Sil2[NumSys][NumConstraint];

        for (int i=0; i<NumSys; i++) {
   	        for (int j=0; j<NumConstraint; j++) {
                Sil2[i][j] = 0;
                num_obs[i][j] = 0;
                R[i][j] = 0;
   			}
   		}

      for (int i=0; i<NumSys; i++) {
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) Z[i][j][d] = 1;
            }
        }

   		for (int i=0; i<NumSys; i++) {

        int demand_index = 0;

   			// generate initial samples
   			double sumY[NumConstraint];
   			double sum_squareY[NumConstraint];

   			for (int j=0; j<NumConstraint; j++) {
   			  sumY[j] = 0;
   				sum_squareY[j] = 0;
   			}

   			for (int n=0; n<Nnot; n++) {
          for (int j=0; j<NumConstraint; j++){
            total_single_obs[j] = 0;
            for (int b=0; b<batch_size; b++) {
              generate_one_obs(i, demand_index);
              total_single_obs[j] += single_sub_obs[j];
            }
            single_obs[j] = total_single_obs[j]/batch_size;
            // printf("obs: %.5f\n",single_obs[1]);
          }

          for (int j=0; j<NumConstraint; j++) {
              sumY[j] += single_obs[j];
              sum_squareY[j] += single_obs[j]*single_obs[j];
              num_obs[i][j] += 1;
            }
          total_obs += 1;
          demand_index += 30;
   			}

            // find continuation region
   			for (int j=0; j<NumConstraint; j++) {
   				Sil2[i][j] = (sum_squareY[j]/(Nnot-1)) - (sumY[j]/(Nnot-1))*(sumY[j]/Nnot);
   				R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
   			}

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];
            for (int j=0; j<NumConstraint; j++) surviveThreshold[j] = NumThreshold;

   			while (surviveConstraint != 0) {

   				for (int j=0; j<NumConstraint; j++) {

   					if (ON[i][j] == 1) {

   						for (int d=0; d<NumThreshold; d++) {
                newq[j][d] = 0;
                  qU[j][d] = q[j][d]*Theta/(q[j][d]*Theta + (1-q[j][d]));
                  qL[j][d] = q[j][d]/(q[j][d] + (1-q[j][d])*Theta);
                  newq[j][d] = (qU[j][d] + qL[j][d])/2;
   							if (ON_l[i][j][d] == 1) {

                  if (((sumY[j]+R[i][j])/num_obs[i][j]) - newq[j][d] <= pow(0.1, 15) ) {
   									Z[i][j][d] = 1;
   									ON_l[i][j][d] = 0;
   									surviveThreshold[j] -= 1;
                  }

                  else if (((sumY[j]-R[i][j])/num_obs[i][j]) - newq[j][d] >= pow(0.1, 15)) {
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

          for (int j=0; j<NumConstraint; j++){
            total_single_obs[j] = 0;
            for (int b=0; b<batch_size; b++) {
              generate_one_obs(i, demand_index);
              total_single_obs[j] += single_sub_obs[j];
            }
            single_obs[j] = total_single_obs[j]/batch_size;
          }

          for (int j=0; j<NumConstraint; j++) {
              sumY[j] += single_obs[j];
              num_obs[i][j] += 1;
          }

          total_obs += 1;
          demand_index += 30;

          for (int j=0; j<NumConstraint; j++) {
              R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
          }
   			}



   			// check whether the decision is correct
   			int cd_for_one_threshold = 1;
            for (int j=0; j<NumConstraint; j++) {
                for (int d=0; d<NumThreshold; d++) {
                  newq[j][d] = 0;
                  qU[j][d] = q[j][d]*Theta/(q[j][d]*Theta + (1-q[j][d]));
                  qL[j][d] = q[j][d]/(q[j][d] + (1-q[j][d])*Theta);
                  newq[j][d] = (qU[j][d] + qL[j][d])/2;
                  // epsilon[j] = (qU[j][d] - qL[j][d])/2;

                	//if (i == 16) printf("%.4f\t%.4f\t%d\n", system_true_value[i][j], q[j][d], Z[i][j][d]);
                  if (system_true_value[i][j] - (newq[j][d] - epsilon[j]) <= pow(0.1, 15)) {
                      if (Z[i][j][d] == 1) cd_for_one_threshold *= 1;
                      else cd_for_one_threshold *= 0;
                  } else if (system_true_value[i][j] - (newq[j][d] + epsilon[j]) >= -pow(0.1, 15)) {
                      if (Z[i][j][d] == 0) cd_for_one_threshold *= 1;
                      else cd_for_one_threshold *= 0;
                  }
                }
            }

   			if (cd_for_one_threshold == 1) {
                final_cd *= 1;
   			} else {
                final_cd *= 0;
            }

   		}

        // for (int i=0; i<NumSys; i++) {
        //     printf("For system %d cosnt1: %d\t%d\t%d\t%d\n", i, Z[i][0][0], Z[i][0][1], Z[i][0][2], Z[i][0][3]);
        //     printf("For system %d const2: %d\t%d\t%d\t%d\n", i, Z[i][1][0], Z[i][1][1], Z[i][1][2], Z[i][1][3]);
        // }

        printf("%.1f\t%.5f\n", total_obs*NumBatch, final_cd);

        write_up();

    }



	return 0;
}
