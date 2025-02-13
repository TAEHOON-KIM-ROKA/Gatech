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
#define NumSys	77
#define NumConstraint	1
#define NumThreshold	4
#define Num_s 20
#define Num_S 20
#define Theta 1.5

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
//double  s10 = 4321, s11 =5432, s12 =6543, s20 =4321, s21 =5432, s22 =6543;
//double  s10 = 11, s11 =12, s12 =13, s20 =21, s21 =22, s22 =23;
//double  s10 = 20202020, s11 =10101010, s12 =30303030, s20 =50505050, s21 =70707070, s22 =90909090;
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
double H[NumConstraint];
double dummies[NumConstraint][NumThreshold];
double qL[NumConstraint][NumThreshold];
double qU[NumConstraint][NumThreshold];
double theta[NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];

int system_value[NumSys][2];
double system_true_value[NumSys][1] = {0};
int true_feasibility[NumSys][NumConstraint][NumThreshold];
double single_obs[NumConstraint];
double demand_list[2000000];

double demand_mean = 25;
double order_cost = 3;
double fixed_order_cost = 32;
double holding_cost = 1;
double penalty_cost = 5;

double total_obs;
double final_cd;

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

  std::ifstream myfile ("real_true_value.txt");
  if (myfile.is_open()) {
    int system_counter = 0;
    while ( system_counter < NumSys && myfile >> prob ) //>> ch >> expected_cost )
    {
      system_true_value[system_counter][0] = prob;
      //system_true_value[system_counter][1] = expected_cost;
      system_counter++;
    }
    myfile.close();
  }

  for (int i = 0; i < NumSys; i++) printf("sys_true: %.10f\n", system_true_value[i][0]);
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

  double total_cost =0;
  double total_fail_prob = 0;
  double LittleS = system_value[system_index][0];
  int BigS = system_value[system_index][1];
  double current_level= BigS, next_level=0, Demand;
  double Cost;

  for(int i=0; i < 12; i++){
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
    }

    current_level = next_level - Demand;
    total_cost += Cost;
  }

  for (int j = 0; j < 1; j++) {
    single_obs[0] = 0;
    //printf("total cost: %.10f\n", total_cost);
    if (total_cost > 1400) {
      single_obs[0] = 1;
    }

    double prn = MRG32k3a();
    for (int d = 0; d < NumThreshold; d++) {
      //double rn = MRG32k3a();
      dummies[j][d] = 0;
      if (prn <= q[j][d]) {
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

  for (int i=0; i<NumSys; i++) {
    for (int j=0; j<NumConstraint; j++) {
      ON[i][j] = 1;

      for (int d=0; d<NumThreshold; d++) {
        ON_l[i][j][d] = 1;
      }
    }
  }

  //epsilon[0] = 0.001;
  //epsilon[1] = 0.5;

// set threshold value
  q[0][0] = 0.01;
  q[0][1] = 0.05;
  q[0][2] = 0.1;
  q[0][3] = 0.2;

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

// struct f {
//     //int idx;
//     //f(int index) : idx(index) {}
//     // 1 threshold. 95%
//     double alpha = 0.05;
//     // >= 2 threshold
//     double k = NumSys;
//     double s = NumConstraint;
//     double beta = (1 - (std::pow( (1-alpha), 1/k))) / (2 * s);

//     double operator()(double ell) {
//         return (1 / (1 + std::pow(theta[0], ell)) - beta);
//     }
// };

// double fsolve() {
//     const boost::uintmax_t maxit = 1000;
//     boost::uintmax_t it = maxit;
//     const double guess = 1.0;
//     const double min = 0.0;
//     const double max = 2.0;
//     boost::math::tools::eps_tolerance<double> tol((std::numeric_limits<double>::digits * 3) / 2);
//     std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(f(), guess, 2.0, false, tol, it);
//     return r.first + (r.second - r.first) / 2;
// };

int main()
{
  read_system_true_value();
  determine_true_feasibility();

  outfile = NULL;
  outfile = fopen("feasibiliy_revision2.out","a");

  double alpha;

  if (NumThreshold == 1){
    alpha = 0.05;
  }
  else {
    alpha = 0.025;
  }

  double k = NumSys;
  double s = NumConstraint;

  printf("alpha: %.4f\n", alpha);
  double beta[NumConstraint];
  for(int j = 0; j < NumConstraint; j++ ){
    //if no CRN
    beta[j] = (1-pow(1-alpha, 1/k)) / (s);   

    //if use CRN
    //beta[j] = (alpha) / (k * s);
  }

  for(int j = 0; j < NumConstraint; j++){
    H[j] = (log ((1/beta[j]) - 1))/(log (Theta));
    H[j] = std::ceil(H[j]);
    printf("H: %.1f\n", H[j]);
  }
  //double eta[NumConstraint];
  //eta[0] = 0.6615;
  //eta[1] = 0.6615;

  for (int l=0; l<NumMacro; l++) {

    configuration();
    total_obs = 0;
    final_cd = 1;

    double num_obs[NumSys][NumConstraint];
    double R[NumSys][NumConstraint];
    double Sil2[NumSys][NumConstraint];

    for (int i=0; i<NumSys; i++) {

      int demand_index = 0;

      // generate initial samples
      double sumY[NumConstraint];
      //double sum_squareY[NumConstraint];
      double sumI[NumConstraint][NumThreshold];

      for (int j=0; j<NumConstraint; j++) {
        sumY[j] = 0;
        for (int d = 0; d < NumThreshold; d++) {
            sumI[j][d] = 0;
        }
      }

      generate_one_obs(i, demand_index);
      demand_index += 30;
      total_obs += 1;

      for (int j = 0; j < NumConstraint; j++) {
        sumY[j] += single_obs[j];
        for (int d = 0; d < NumThreshold; d++) {
          sumI[j][d] += dummies[j][d];
        }
        num_obs[i][j] += 1;     //각 constraint마다 존재하는 obs의 갯수
      }

      int surviveConstraint = NumConstraint;
      int surviveThreshold[NumConstraint];
      for (int j=0; j<NumConstraint; j++) surviveThreshold[j] = NumThreshold;

      while (surviveConstraint != 0) {

        for (int j=0; j<NumConstraint; j++) {

          if (ON[i][j] == 1) {

            for (int d=0; d<NumThreshold; d++) {

              if (ON_l[i][j][d] == 1) {

                if (sumY[j] - sumI[j][d] <= -H[j]) {
                  Z[i][j][d] = 1;
                  ON_l[i][j][d] = 0;
                  surviveThreshold[j] -= 1;
                }

                else if (sumY[j] - sumI[j][d] >= H[j]) {
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

        generate_one_obs(i, demand_index);
        demand_index += 30;
        total_obs += 1;

        for (int j = 0; j < NumConstraint; j++) {
          sumY[j] += single_obs[j];
          for (int d = 0; d < NumThreshold; d++) {
            sumI[j][d] += dummies[j][d];
          }
          num_obs[i][j] += 1;
        }
        //printf("sumY: %.10f\n", sumY[0]);
        //printf("sumI: %.10f\n", sumI[0][0]);
      }

      // check whether the decision is correct
      int cd_for_one_threshold = 1;
        for (int j=0; j<NumConstraint; j++) {
            
          for (int d=0; d<NumThreshold; d++) {

            if (system_true_value[i][j] <= (q[j][d] / (q[j][d] + (1 - q[j][d]) * Theta))) {
              if (Z[i][j][d] == 1) cd_for_one_threshold *= 1;
              else cd_for_one_threshold *= 0;
            } 
            else if (system_true_value[i][j] >= (q[j][d] * Theta / ((1 - q[j][d]) + q[j][d] * Theta))) {
              if (Z[i][j][d] == 0) cd_for_one_threshold *= 1;
              else cd_for_one_threshold *= 0;
            }
            else {
              cd_for_one_threshold *= 1;
            }
          }
        }

      if (cd_for_one_threshold == 1) {
        final_cd *= 1;
      } else {
        final_cd *= 0;
      }

    }

    // for (int i=0; i<NumSys; i++){
    //   for (int j=0; j<NumConstraint; j++){
    //     for (int d=0; d<NumThreshold; d++){
    //       printf("%d\t%d\t%d\t%d\n", i, j, d, Z[i][j][d]);
    //     }
    //   }
    // }

    printf("%.1f\t%.5f\n", total_obs, final_cd);

    write_up();

  }
	return 0;
}