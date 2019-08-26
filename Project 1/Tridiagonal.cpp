#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
using namespace std;

ofstream ofile;

//Function for f and exact solution, respectively
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double u_analytic(double x){return (exp(-10.0)-1.0)*x+1.0-exp(-10.0*x);}

int main(int argc, char *argv[]){
  int max_m;            //maximum value of m (n = 10^m)
  string fname;         //beginning of filename; filenames will be "fname_m"
  int solution_method;  //1 for general method, 2 for special case
  if(argc <= 3){
    cout << "No filename, m and/or solution method; read filename, max. " "value of m (n=10^m) and solution method (1 for general case, 2 for "
    "special) on the same line." << endl;
    exit(1);
  }
  else if (atoi(argv[3]) != 1 and atoi(argv[3]) != 2){
    cout << "Invalid value for solution method; input 1 for general case, "
    "2 for special." << endl;
    exit(1);
  }
  else{
    fname = argv[1];
    max_m = atoi(argv[2]);
    solution_method = atoi(argv[3]);
  }
  double *t = new double[max_m];           //Time per n
  double *max_rel_err = new double[max_m]; //Maximum rel. error per n
  for (int m = 1; m <= max_m; m++){
    int n = int(pow(10.,m));
    string outfile = fname;
    outfile.append("_"+to_string(m));

    //Setup of variables. All arrays are set to length n+2 for convenience;
    //for a, b, and c elements 0 and n+1 are not used (element n is not used
    //either for a and c)
    clock_t start, finish;
    start = clock();

    double h = 1.0/(n+1);                  //step size
    double h2 = h*h;
    double *a = new double[n+2];       //"below" diagonal
    double *b = new double[n+2];       //diagonal
    double *c = new double[n+2];       //"above" diagonal
    double *x = new double[n+2];
    double *g = new double[n+2];       //numerical solution, initially equal
                                       //to h^2*f(x_i) for i=1,2,...,n
    double *u = new double[n+2];       //exact solution
    for (int i = 0; i <= n+1; i++){
      x[i] = i*h;
      u[i] = u_analytic(x[i]);
      a[i] = c[i] = -1.;
      b[i] = 2.;
      g[i] = h2*f(x[i]);
    }
    g[0] = g[n+1] = 0.;         //Boundary conditions on u(x)
    u[0] = u[n+1] = 0.;         //Removes round-off error in analytic expression
    //Solving the equation (with either method):
    if (solution_method == 1){
      for (int i = 2; i <= n; i++){
        double a_b = a[i-1]/b[i-1];
        b[i] -= a_b*c[i-1];
        g[i] -= a_b*g[i-1];
      }
      g[n] /= b[n];
      for (int i=n-1; i >= 1; i--){
        g[i] = (g[i]-c[i]*g[i+1])/b[i];  //solution
      }
    }
    else if (solution_method == 2){
      cout << "Under construction!" << endl;
      exit(1);
    }
    finish = clock();
    t[m-1] = (finish-start)/double(CLOCKS_PER_SEC);
    double *rel_err = new double[n+2];
    //Writing solution to file (for m up to 3):
    if(m <= 3){
      ofile.open(outfile);
      ofile << setiosflags(ios::showpoint);
      ofile << "             x:   Numerical solution:  "
      "Exact solution:  Absolute error:  Relative error:" << endl;
      for (int i = 0; i <= n+1; i++){
        double err = g[i]-u[i];
        if(err == 0.) rel_err[i] = 0.;
        else rel_err[i] = fabs(err)/u[i];
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(22) << setprecision(8) << g[i];
        ofile << setw(17) << setprecision(8) << u[i];
        ofile << setw(17) << setprecision(8) << err;
        ofile << setw(17) << setprecision(8) << rel_err[i] << endl;
      }
      ofile.close();
    }
    else{
      for (int i = 0; i <= n+1; i++){
        double err = g[i]-u[i];
        if(err == 0.) rel_err[i] = 0.;
        else rel_err[i] = fabs(err)/u[i];
    }}
    //Find maximum relative error for this n:
    max_rel_err[m-1] = *max_element(rel_err, rel_err + n+2);

    delete [] a; delete [] b; delete [] c;
    delete [] x; delete [] g; delete [] u;
    }
  //Writing times to file:
  ofile.open(fname+"_time");
  ofile << setiosflags(ios::showpoint);
  ofile << "  n:           t [s]:" << endl;
  for (int m=1; m <=  max_m; m++){
    ofile << "10^" << m;
    ofile << setw(16) << setprecision(8) << t[m-1] << endl;
  }
  ofile.close();
  //Writing relative errors to file:
  ofile.open(fname+"_error");
  ofile << setiosflags(ios::showpoint);
  ofile << "        n:   max. rel. error:" << endl;
  for (int m=1; m <=  max_m; m++){
    ofile << setw(10) << int(pow(10,m));
    ofile << setw(16) << setprecision(8) << max_rel_err[m-1] << endl;
  }
  ofile.close();

  return 0;
}
