#include <iostream>
#include <vector>
#include <cmath>
#include <utility> // for std::pair
#include <string>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <random>
#include <list>

using namespace std;

// user inputs for N_0, number of macro-rep, number of systems, number of constraints
// and number of thresholds of all constraint (if constraints have different number
// of threshods, then input the maximum number of threshods and adjust the actual
// number of thresholds each constraint later in the code)
#define NumMacro 10000
#define NumSys	77
#define NumConstraint	1
#define NumThreshold	4
#define Num_s 20
#define Num_S 20
#define Theta 1.2

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
double demand_list[15000000];

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

  for (int i=0; i<15000000; i++) {
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
    //Demand = poisson(demand_mean);
    Demand = demand_list[demand_index+i];

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


// 상관계수를 계산하는 함수
double calculateCorrelation(const vector<double>& x, const vector<double>& y) {
    size_t n = x.size();
    double sumX = accumulate(x.begin(), x.end(), 0.0);
    double sumY = accumulate(y.begin(), y.end(), 0.0);
    double meanX = sumX / n;
    double meanY = sumY / n;
    
    double sumXY = 0, sumX2 = 0, sumY2 = 0;
    for (size_t i = 0; i < n; ++i) {
        sumXY += (x[i] - meanX) * (y[i] - meanY);
        sumX2 += pow(x[i] - meanX, 2);
        sumY2 += pow(y[i] - meanY, 2);
    }
    
    return sumXY / (sqrt(sumX2) * sqrt(sumY2));
}

void generateSystemData(vector<vector<double>>& total_costs, vector<vector<double>>& single_obs_values, int numSystems, int observationsPerSystem) {
    // total_costs와 single_obs_values 벡터를 초기화합니다.
    total_costs.resize(numSystems, vector<double>(observationsPerSystem, 0));
    single_obs_values.resize(numSystems, vector<double>(observationsPerSystem, 0));
    
    // 각 시스템별로 total_cost와 single_obs 값을 생성합니다.
    for (int systemIndex = 0; systemIndex < numSystems; ++systemIndex) {
        int demand_index = 0;
        for (int obsIndex = 0; obsIndex < observationsPerSystem; ++obsIndex) {
            // 여기에서 각 시스템별로 total_cost와 single_obs 값을 생성하는 로직을 구현합니다.
            // 예시: total_costs[systemIndex][obsIndex] = (생성 로직);
            // 예시: single_obs_values[systemIndex][obsIndex] = (생성 로직);
            // 주어진 코드의 `generate_one_obs` 함수 로직 참조
            double total_cost =0;
            double total_fail_prob = 0;
            double LittleS = system_value[systemIndex][0];
            int BigS = system_value[systemIndex][1];
            double current_level= BigS, next_level=0, Demand;
            double Cost;

            for(int i=0; i < 12; i++){
                Cost = 0;
                // Demand = poisson(demand_mean);
                Demand = demand_list[demand_index+i];

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

            // for (int j = 0; j < 1; j++) {
            //     single_obs[0] = 0;
            //     //printf("total cost: %.10f\n", total_cost);
            //     if (total_cost > 1400) {
            //     single_obs[0] = 1;
            //     }
            // }
            // printf("System %.1d total cost: %.1f\n", systemIndex, total_cost);
            total_costs[systemIndex][obsIndex] = total_cost;
            // single_obs_values[systemIndex][obsIndex] = single_obs[0];
            demand_index += 12;
        }
        printf("System %.1d obs generated completely.\n", systemIndex);
    }
}

int main() {

    generate_demand();
    configuration();
    const int ObservationsPerSystem = 1000000; // 각 시스템별 관측값의 수

    vector<vector<double>> total_costs;
    vector<vector<double>> single_obs_values;
    
    // 시스템별 데이터 생성
    generateSystemData(total_costs, single_obs_values, NumSys, ObservationsPerSystem);
    
    // 각 시스템 쌍에 대해 total_cost의 상관계수 계산
    cout << "Correlation of total_costs between systems:" << endl;
    for (int i = 0; i < NumSys; ++i) {
        for (int j = i + 1; j < NumSys; ++j) {
            double correlation = calculateCorrelation(total_costs[i], total_costs[j]);
            cout << i << " " << j << " " << correlation << endl;
        }
    }
    
    // 각 시스템 쌍에 대해 single_obs의 상관계수 계산
    // cout << "\nCorrelation of single_obs between systems:" << endl;
    // for (int i = 0; i < NumSys; ++i) {
    //     for (int j = i + 1; j < NumSys; ++j) {
    //         double correlation = calculateCorrelation(single_obs_values[i], single_obs_values[j]);
    //         cout << i << " " << j << " " << correlation << endl;
    //     }
    // }
    
    return 0;
}