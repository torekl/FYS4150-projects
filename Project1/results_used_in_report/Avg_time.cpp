#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <numeric>
#include <armadillo>
using namespace std;
using namespace arma;
ofstream ofile;

double avg(double t[],int n){
  double s = 0.;
  return accumulate(t, t+n, s)/double(n);
}


int main(int argc, char *argv[]){
  mat t_general(10,7);
  mat t_special(10,7);
  mat t_LU(10,4);

  //Times from general algorithm:
  for (int run = 1; run <= 10; run++){
    ifstream infile("general_run"+to_string(run)+"_time");
    string headline;
    getline(infile,headline); //skip first line
    string n;
    double t;
    int m = 0;
    while(infile >> n >> t){
      t_general(run-1,m) = t;
      m++;
    }
    infile.close();
  }
  ofile.open("general_avg_time");
  ofile << setiosflags(ios::showpoint);
  ofile << "  n:           t [s]:" << endl;
  for (int m=1; m <=  7; m++){
    double *row_m = t_general.colptr(m-1); //Times for 10^m integration points
    double t_avg = avg(row_m,10);
    ofile << "10^" << m;
    ofile << setw(16) << setprecision(6) << t_avg << endl;
  }
  ofile.close();

  //Times from special algorithm:
  for (int run = 1; run <= 10; run++){
    ifstream infile("special_run"+to_string(run)+"_time");
    string headline;
    getline(infile,headline); //skip first line
    string n;
    double t;
    int m = 0;
    while(infile >> n >> t){
      t_special(run-1,m) = t;
      m++;
    }
    infile.close();
  }
  ofile.open("special_avg_time");
  ofile << setiosflags(ios::showpoint);
  ofile << "  n:           t [s]:" << endl;
  for (int m=1; m <=  7; m++){
    double *row_m = t_special.colptr(m-1); //Times for 10^m integration points
    double t_avg = avg(row_m,10);
    ofile << "10^" << m;
    ofile << setw(16) << setprecision(6) << t_avg << endl;
  }
  ofile.close();

  //Times from LU decomposition:
  for (int run = 1; run <= 10; run++){
    ifstream infile("LU_run"+to_string(run)+"_time");
    string headline;
    getline(infile,headline); //skip first line
    string n;
    double t;
    int m = 0;
    while(infile >> n >> t){
      t_LU(run-1,m) = t;
      m++;
    }
    infile.close();
  }
  ofile.open("LU_avg_time");
  ofile << setiosflags(ios::showpoint);
  ofile << "  n:           t [s]:" << endl;
  for (int m=1; m <=  4; m++){
    double *row_m = t_LU.colptr(m-1); //Times for 10^m integration points
    double t_avg = avg(row_m,10);
    ofile << "10^" << m;
    ofile << setw(16) << setprecision(6) << t_avg << endl;
  }
  ofile.close();
  return 0;
}
