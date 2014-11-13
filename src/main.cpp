#include <iostream>
#include <cmath>
#include <omp.h>
#include "SOR.h"
using namespace std;


double righthandside(double x,double y){
  //x = 0; y = 0
  return 0;
}

double leftBC(const double& x, const double& y){
  return 0;
}

double rightBC(const double& x, const double& y){
  return 0;
}

double topBC(const double& x, const double& y){
  return sin(M_PI*x)*exp(-M_PI);
}

double bottomBC(const double& x, const double& y){
  return sin(M_PI*x);
}

double exactsol(const double& x, const double& y){
  return sin(M_PI*x)*exp(-M_PI*y);
}



int main()
{
  int nthreads = 1;
  //cout << "Enter number of threads : ";
  //cin >> nthreads;
  omp_set_num_threads(nthreads);

  PoissonEq pe(1,1,righthandside,0,1,0,1);
  Grid grid(100,100,pe);
  BoundaryConditions bc;
  // set boundary conditions
  bc.SetLeftBC(leftBC);
  bc.SetRightBC(rightBC);
  bc.SetTopBC(topBC);
  bc.SetBottomBC(bottomBC);

  SOR solver(&pe,&bc,&grid);
  solver.setMaxIterations(1e8);
  solver.setRelaxationParameter(1.9);
  //solver.set_exactsol_fun(exactsol);
  //solver.exact_solution();

  double tstart = omp_get_wtime();
  solver.solve_sor_omp();
  double tend = omp_get_wtime();

  solver.WriteSolution("sol.dat");

  cout << "Time taken = " << tend - tstart << endl;
  return 0;
}
