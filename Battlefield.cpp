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

#define NumMacro 1
#define NumSys	84
#define NumConstraint	2
#define NumThreshold	5
#define Nnot 20
#define Theta   1.2
#define Epsilon 2.5

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
// double  s10 = 43, s11 =54, s12 =65, s20 =43, s21 =54, s22 =65;
// double  s10 = 4321111, s11 =1115432, s12 =1116543, s20 =4321111, s21 =1115432, s22 =6543111;
//double  s10 = 43221, s11 =54332, s12 =65443, s20 =43321, s21 =54532, s22 =61543;
//double  s10 = 1010, s11 =10, s12 =101, s20 =2001, s21 = 202, s22 = 202;

//double normal(double rmean, double rvar);
//double generate_multiNormal(int numConstraint, int case_index);
double configuration(void);
//double generate_Bernoulli(int numConstraint, int case_index);
//int read_chol_matrix(void);
//float correlationCoefficient(std::vector<double> X, std::vector<double> Y, int n);

double mean_value[NumSys][NumConstraint];
//double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys];

double observations[NumSys][NumConstraint];
double H[NumSys][NumConstraint];
double dummies[NumThreshold][NumConstraint];
double q[NumThreshold][NumConstraint];
double qL[NumThreshold][NumConstraint];
double qU[NumThreshold][NumConstraint];
double newq[NumThreshold][NumConstraint];
double alpha;
double k;
double s;
double beta;
double eta[NumConstraint];
double theta[NumConstraint];
double thetas[NumThreshold][NumConstraint];
double epsilon[NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];
double system_true_value[NumSys][2] = {0};

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

int read_system_true_value(void) {
    double prob;
    double expected_cost;
    char comma;

    std::ifstream myfile("RealTrueMean.csv");
    if (myfile.is_open()) {
        int system_counter = 0;
        while (system_counter < NumSys && myfile >> prob >> comma >> expected_cost) {
            system_true_value[system_counter][0] = prob;
            system_true_value[system_counter][1] = expected_cost;
            system_counter++;
        }
        myfile.close();
    }

    // for (int i = 0; i < NumSys; i++) {
    //     printf("sys_true_prob: %.4f\n", system_true_value[i][0]);
    //     printf("sys_true_cost: %.4f\n", system_true_value[i][1]);
    // }
    return 0;
}

struct Observation {
    double observation1;
    double observation2;
};

class CSVReader {
public:
    CSVReader(const string& filename) : file(filename), isEOF(false) {}

    Observation readNextObservation() {
        Observation obs;

        if (!file.is_open()) {
            cerr << "Failed to open file" << endl;
            return obs;  // Return empty observations
        }

        string line;
        if (getline(file, line)) {
            istringstream line_stream(line);
            string cell;
            int cell_num = 0;
            while (getline(line_stream, cell, ',')) {
                if (cell_num == 7) {
                    obs.observation1 = stod(cell);
                } else if (cell_num == 8) {
                    obs.observation2 = stod(cell);
                }
                cell_num++;
            }
        } else {
            isEOF = true;  // Mark EOF
        }

        return obs;
    }

    bool isEOFReached() const {
        return isEOF;
    }

private:
    ifstream file;
    bool isEOF;
};

double configuration(void) {

    for (int i = 0; i < NumSys; i++) {
        for (int j = 0; j < NumConstraint; j++) {
            ON[i][j] = 1;
            for (int d = 0; d < NumThreshold; d++) {
                ON_l[i][j][d] = 1;
            }
        }
    }

    /*for (int j = 0; j < NumConstraint; j++) {
        epsilon[0] = 0.05;
        epsilon[1] = 0.05;
    }*/

    // two constraints and two thresholds
    q[0][0] = 0.05;
    q[1][0] = 0.1;
    q[2][0] = 0.15;
    q[3][0] = 0.2;
    q[4][0] = 0.25;
    q[0][1] = 620;
    q[1][1] = 640;
    q[2][1] = 660;
    q[3][1] = 680;
    q[4][1] = 700;

    // q[0][0] = 0.01;
    // q[1][0] = 0.05;
    // q[2][0] = 0.1;
    // q[3][0] = 0.15;
    // q[4][0] = 0.2;
    // q[5][0] = 0.3;
    // q[0][1] = 600;
    // q[1][1] = 620;
    // q[2][1] = 640;
    // q[3][1] = 660;
    // q[4][1] = 680;
    // q[5][1] = 700;

    theta[0] = Theta;
    epsilon[1] = Epsilon;

    return 0;
}


struct f {
    //int idx;
    //f(int index) : idx(index) {}
    // 1 threshold. 95%
    //double alpha = 0.05;
    // >= 2 threshold
    double alpha = 0.05;
    double k = NumSys;
    double s = NumConstraint;
    double beta = (1 - (std::pow( (1-alpha), 1/k))) / (2*s);

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

int main()
{

    read_system_true_value();

    double total_obs;
    double final_cd;   
    double correct_decision;   
    double overall_obs; 

    double n0_counter;

    correct_decision = 0;
    overall_obs = 0;

    for (int l = 0; l < NumMacro; l++) {   

        configuration();   
        total_obs = 0;      
        final_cd = 1;  
        n0_counter = 0;    

        double num_obs[NumSys][NumConstraint];
        double R[NumSys][NumConstraint];
        double Sil2[NumSys][NumConstraint];

        for (int i = 0; i < NumSys; i++) {
            for (int j = 0; j < NumConstraint; j++) {
                //if (j = 0){
                    H[i][j] = std::ceil(fsolve());
                    //printf("theta: %.10f\n", theta[j]);
                    //printf("H: %.10f\n", H[i][j]);

                    num_obs[i][0] = 0;
                //}
                if (j = 1){
                    double alpha = 0.05;
                    double k = NumSys;
                    double s = NumConstraint;
                    double n0 = Nnot;

                    double beta = (1 - (std::pow( (1-alpha), 1/k))) / (2*s);
                    //printf("beta: %.10f\n", beta);
                    eta[1] = (pow(2 * beta, -1 * (2 / (n0 - 1))) - 1) / 2;
                    //printf("eta: %.10f\n", eta[1]);

                    Sil2[i][1] = 0;
                    num_obs[i][1] = 0;
                }
            }
        }

        for (int i = 0; i < NumSys; i++) {

            n0_counter = 0;

            string filename = "data/s" + to_string(i + 1) + ".csv";
            CSVReader reader(filename);

            Observation obs = reader.readNextObservation();
            observations[i][0] = obs.observation1;
            observations[i][1] = obs.observation2;

            //printf("obs value 1: %.4f\n", observations[i][0]);
            //printf("obs value 2: %.4f\n", observations[i][1]);

            n0_counter += 1;

            // generate initial samples
            double sumY[NumConstraint];  
            double sum_squareY[NumConstraint];  
            double sumI[NumThreshold][NumConstraint];   

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

            for (int j = 0; j < NumConstraint; j++) {
                sumY[j] = 0;
                sum_squareY[j] = 0;
                
            }

            for (int d = 0; d < NumThreshold; d++) {
                    sumI[d][0] = 0;
            }

            //dummy generator
            double rn = MRG32k3a();
            for (int d = 0; d < NumThreshold; d++) {
                dummies[d][0] = 0;
                if (rn <= q[d][0]) {
                    dummies[d][0] = 1;
                }
            }

            total_obs += 1; 

            for (int j = 0; j < NumConstraint; j++) {
                
                sumY[j] += observations[i][j];
                sum_squareY[j] += observations[i][j] * observations[i][j];
                num_obs[i][j] += 1;
            }

            for (int j = 0; j < NumConstraint; j++) surviveThreshold[j] = NumThreshold;

            for (int d = 0; d < NumThreshold; d++) {
                    sumI[d][0] += dummies[d][0];
            }

            while (surviveConstraint != 0) {

                //for (int j = 0; j < NumConstraint; j++) {

                    //if (j = 0){
                        int j = 0;

                        if (ON[i][j] == 1) {

                            for (int d = 0; d < NumThreshold; d++) { 

                                if (ON_l[i][j][d] == 1) {

                                    if (sumY[j] - sumI[d][j] <= -H[i][j]) { 
                                        Z[i][j][d] = 1;     
                                        ON_l[i][j][d] = 0; 
                                        surviveThreshold[j] -= 1;
                                    }

                                    else if (sumY[j] - sumI[d][j] >= H[i][j]) { 
                                        Z[i][j][d] = 0;    
                                        ON_l[i][j][d] = 0;    
                                        surviveThreshold[j] -= 1;
                                    }
                                }   
                            }

                            if (surviveThreshold[j] == 0) { 
                                ON[i][j] = 0;
                                surviveConstraint -= 1;
                                printf("s %d BeRF end with %.0f data \n", i+1, num_obs[i][0]);
                            }
                        }
                    //}

                    //else if (j = 1){
                        j = 1;

                        if (n0_counter >= Nnot){

                            if (ON[i][j] == 1) {    //만약 그 시스템에서 그 constraint가 아직 검사가 안됐다면

                                for (int d = 0; d < NumThreshold; d++) {    //threshold까지 내려가서 검사

                                    if (ON_l[i][j][d] == 1) {   //검사 안됐었으면

                                        if ((sumY[j] + R[i][j]) / num_obs[i][j] <= q[d][j]) {     // feasible 조건
                                            Z[i][j][d] = 1;     // i 시스템의 j 제약의 d번째 threshold에 대해 feasible
                                            ON_l[i][j][d] = 0;  // 해당 threshold를 검사한거로 변경
                                            surviveThreshold[j] -= 1;
                                        }

                                        else if ((sumY[j] - R[i][j]) / num_obs[i][j] >= q[d][j]) {   // infeasible 조건
                                            Z[i][j][d] = 0;     // i 시스템의 j 제약의 d번째 threshold에 대해 infeasible
                                            ON_l[i][j][d] = 0;      // 해당 threshold를 검사한거로 변경
                                            surviveThreshold[j] -= 1;
                                        }
                                    }
                                }

                                if (surviveThreshold[j] == 0) {     //더이상 검사안한 threshold가 없으면 해당 constraint도 검사한거로 변경
                                    ON[i][j] = 0;
                                    surviveConstraint -= 1;
                                    printf("s %d RF end with %.0f data \n", i+1, num_obs[i][1]);
                                }
                            }
                        }
                    //}
                //}

                if (surviveConstraint == 0) break; 

                else if (reader.isEOFReached()) { 
                    
                    printf("****system %d needs more data.**** \n", i+1);

                    break;
                }

                Observation obs = reader.readNextObservation();
                observations[i][0] = obs.observation1;
                observations[i][1] = obs.observation2;

                //printf("obs value 1: %.4f\n", observations[i][0]);
                //printf("obs value 2: %.4f\n", observations[i][1]);

                //dummy generator
                double rn = MRG32k3a();
                for (int d = 0; d < NumThreshold; d++) {
                    dummies[d][0] = 0;
                    if (rn <= q[d][0]) {
                        dummies[d][0] = 1;
                    }
                }
                
                n0_counter += 1;
                total_obs += 1; 

                for (int j = 0; j < NumConstraint; j++) {
                    sumY[j] += observations[i][j];
                    sum_squareY[j] += observations[i][j] * observations[i][j];
                    num_obs[i][j] += 1;
                }

                for (int d = 0; d < NumThreshold; d++) {
                        sumI[d][0] += dummies[d][0];
                }

                if (n0_counter == Nnot){
                    Sil2[i][1] = (sum_squareY[1] / (Nnot - 1)) -(sumY[1] / (Nnot - 1)) * (sumY[1] / Nnot);
                }

                R[i][1] = maxfn(0, (Nnot - 1) * Sil2[i][1] * (eta[1]) / epsilon[1] - epsilon[1] * num_obs[i][1] / 2);
                //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];

            }

            // check whether the decision is correct
            int cd_for_one_threshold = 1;   
            for (int d = 0; d < NumThreshold; d++) {  
                //for (int j = 0; j < NumConstraint; j++) { 
                    //if (j = 0) {
                        //mean_value[i][j] = sumY[j] / total_obs;
                        int j = 0;
                        if (system_true_value[i][j] <= (q[d][j] / (q[d][j] + (1 - q[d][j]) * theta[j]))) {  
                            if (Z[i][j][d] == 1) {     
                                cd_for_one_threshold *= 1; 
                                //printf("s %d cons1 Correct\n", i+1);
                            }
                            else {
                                cd_for_one_threshold *= 0;
                                //printf("s %d cons1 Incorrect\n", i+1);  
                            }
                        }
                        else if (system_true_value[i][j] >= (q[d][j] * theta[j] / ((1 - q[d][j]) + q[d][j] * theta[j]))) {   
                            if (Z[i][j][d] == 0) {  
                                cd_for_one_threshold *= 1;
                                //printf("s %d cons1 Correct\n", i+1);  
                            }
                            else {
                                cd_for_one_threshold *= 0;
                                //printf("s %d cons1 Incorrect\n", i+1); 
                            }
                        }
                        else {
                            cd_for_one_threshold *= 1;
                            //printf("s %d cons1 Correct\n", i+1); 
                        }
                    //}
                    //else if (j = 1) {
                        //mean_value[i][j] = sumY[j] / num_obs[i][j];
                        j = 1;
                        if (system_true_value[i][j] <= q[d][j] - epsilon[j]) {     //tolerance level을 뺀 q보다 mean value가 작다면,
                            if (Z[i][j][d] == 1) {      // feasible하게 판단했다면,
                                cd_for_one_threshold *= 1;  // correct decision
                                //printf("s %d cons2 Correct\n", i+1);
                            }
                            else {
                                cd_for_one_threshold *= 0;  // incorrect decision
                                //printf("s %d cons2 Incorrect\n", i+1);
                            }
                        }
                        else if (system_true_value[i][j] >= q[d][j] + epsilon[j]) {    //tolerance level을 더한 q보다 mean value가 크다면,
                            if (Z[i][j][d] == 0) {  //만약 infeasible하게 판단했다면, 
                                cd_for_one_threshold *= 1;  //correct decision
                                //printf("s %d cons2 Correct\n", i+1);
                            }
                            else {
                                cd_for_one_threshold *= 0;  //incorrect
                                //printf("s %d cons2 Incorrect\n", i+1);
                            }
                        }
                        else {
                            cd_for_one_threshold *= 1;  //Acceptable하니깐 그냥 correct decision이라고 하자.
                            //printf("s %d cons2 Correct\n", i+1);
                        }
                    //}
                //}
            }

            if (cd_for_one_threshold == 1) {   
                final_cd *= 1;  
            }
            else {
                final_cd *= 0;  
            }

        }

        if (final_cd == 1) {    //최종 1이면 correct decision의 개수 하나 올려주기
            correct_decision++;
        }

        overall_obs += total_obs;   //각 iter마다 했던 obs 갯수를 다 더함

    }

    // for (int i=0; i<NumSys; i++){
    //     for (int j=0; j<NumConstraint; j++){
    //         for (int d=0; d<NumThreshold; d++){
    //             printf("%d\t%d\t%d\t%d\n", i+1, j, d, Z[i][j][d]);
    //         }
    //     }
    // }

    printf("Overall PCD: %.10f\n", correct_decision / NumMacro);
    printf("Overall OBS: %.4f\n", overall_obs / NumMacro);

    return 0;
}
