#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
using namespace std;
using namespace std::chrono;
using namespace arma;

ofstream ofile;

//Function for f and exact solution, respectively
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double u_analytic(double x){return (exp(-10.0)-1.0)*x+1.0-exp(-10.0*x);}

int main(int argc, char *argv[]){
  int max_m;            //maximum value of m (n = 10^m)
  string fname;         //beginning of filename; filenames will be "fname_m"
  int solution_method;  //1 for general method, 2 for special case
  if(argc <= 3){
    cout << "No filename, m and/or solution method; read filename, max. "
    "value of m (n=10^m) and solution method (1 for general tridiagonal "
    "case, 2 for special, 3 for LU decomp.) on the same line." << endl;
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

    //Setup of variables. All arrays are set to length n+2 for convenience,
    //(so that the relevant elements have indices 1-n)even though for some
    //the first/last elements are not used
    time_point<high_resolution_clock> start, finish;
    double h = 1.0/(n+1);                  //step size
    double h2 = h*h;
    double *a = new double[n+2];       //"below" diagonal
    double *b = new double[n+2];       //diagonal
    double *c = new double[n+2];       //"above" diagonal
    double *x = new double[n+2];
    double *g = new double[n+2];       //h^2*f(x_i) for i=1,2,...,n
    double *u_exact = new double[n+2]; //exact solution
    double *b_new = new double[n+2];
    double *g_new = new double[n+2];
    double *u_num = new double[n+2];
    for (int i = 0; i <= n+1; i++){
      x[i] = i*h;
      u_exact[i] = u_analytic(x[i]);
      a[i] = c[i] = -1.;
      b[i] = 2.;
      g[i] = h2*f(x[i]);
    }
    u_num[0] = u_num[n+1] = 0.;
    u_exact[0] = u_exact[n+1] = 0.;    //Removes round-off error
                                       //in analytic expression
    //Solving the equation (with either method):

    //General tridiagonal solver:
    if (solution_method == 1){
      b_new[1] = b[1]; g_new[1] = g[1];
      start = high_resolution_clock::now();
      for (int i = 2; i <= n; i++){
        double a_b = a[i-1]/b_new[i-1];
        b_new[i] = b[i] - a_b*c[i-1];
        g_new[i] = g[i] - a_b*g_new[i-1];
      }
      u_num[n] = g_new[n]/b_new[n];
      for (int i=n-1; i >= 1; i--){
        u_num[i] = (g_new[i]-c[i]*u_num[i+1])/b_new[i];  //solution
      }
    }

    //Specialized tridiagonal solver:
    else if (solution_method == 2){
      for (int i = 1; i <= n; i++) b_new[i] = (i+1.)/i;
      g_new[1] = g[1];
      start = high_resolution_clock::now();
      for (int i = 2; i <= n; i++){
        g_new[i] = g[i] + g_new[i-1]/b_new[i-1];
      }
      u_num[n] = g_new[n]/b_new[n];
      for (int i=n-1; i >= 1; i--){
        u_num[i] = (g_new[i]+u_num[i+1])/b_new[i];  //solution
      }
    }

    //LU decomposition:
    else if (solution_method == 3){
      mat A(n,n);
      mat L, U, P;
      vec g_vec(n);
      //setup of A and g:
      for (int i = 0; i <= n-1; i++){
        for (int j = 0; j <= n-1; j++){
          if (i == j) A(i,j) = 2.;
          else if (i == j + 1 or i == j - 1) A(i,j) = -1.;
          else A(i,j) = 0;
        }
        g_vec[i] = g[i+1];
      }
      start = high_resolution_clock::now();
      lu(L, U, P, A);
      vec y = solve(L,g_vec);
      vec solution = solve(U,y);
      for (int i = 1; i <= n; i++){
        u_num[i] = solution[i-1];
      }
    }
    else{
      cout << "Invalid value for solution method; input 1 for general "
      "tridiagonal case, 2 for special, 3 for LU decomp." << endl;
      exit(1);
    }
    finish = high_resolution_clock::now();
    duration<double> elapsed = finish-start;
    t[m-1] = elapsed.count();
    double *rel_err = new double[n+2];
    //Writing solution to file (for m up to 3):
    if(m <= 3){
      ofile.open(outfile);
      ofile << setiosflags(ios::showpoint);
      ofile << "             x:   Numerical solution:  "
      "Exact solution:  Absolute error:  Relative error:" << endl;
      for (int i = 0; i <= n+1; i++){
        double err = u_num[i]-u_exact[i];
        if(err == 0.) rel_err[i] = 0.;
        else rel_err[i] = fabs(err)/u_exact[i];
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(22) << setprecision(8) << u_num[i];
        ofile << setw(17) << setprecision(8) << u_exact[i];
        ofile << setw(17) << setprecision(8) << err;
        ofile << setw(17) << setprecision(8) << rel_err[i] << endl;
      }
      ofile.close();
    }
    else{
      for (int i = 0; i <= n+1; i++){
        double err = u_num[i]-u_exact[i];
        if(err == 0.) rel_err[i] = 0.;
        else rel_err[i] = fabs(err)/u_exact[i];
    }}
    //Find maximum relative error for this n:
    max_rel_err[m-1] = *max_element(rel_err, rel_err + n+2);

    delete [] a; delete [] b; delete [] c;
    delete [] x; delete [] g; delete [] u_exact;
    delete [] b_new; delete [] g_new; delete [] u_num;
    }
  //Writing times to file:
  ofile.open(fname+"_time");
  ofile << setiosflags(ios::showpoint);
  ofile << "  n:           t [s]:" << endl;
  for (int m=1; m <=  max_m; m++){
    ofile << "10^" << m;
    ofile << setw(16) << setprecision(6) << t[m-1] << endl;
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
