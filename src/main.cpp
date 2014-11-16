#include <iostream>
#include <cmath>
#include <omp.h>
#include "mpi.h"
#include "sor.hpp"
using namespace std;

/* Error Check */
#define CHKERRQ(n) if(n != MPI_SUCCESS) printf("Error!! Check line number : %d\n",__LINE__)

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


int main(int argc, char **argv)
{
  int ierr;

  ierr = MPI_Init(&argc,&argv); CHKERRQ(ierr);


  PoissonEq pe(1,1,righthandside,0,1,0,1);
  Grid grid;
  grid.SetGridSize(100,100);
  grid.SetEquation(pe);
  BoundaryConditions bc;
  // set boundary conditions
  bc.SetLeftBC(leftBC);
  bc.SetRightBC(rightBC);
  bc.SetTopBC(topBC);
  bc.SetBottomBC(bottomBC);

  SOR solver(&pe,&bc,&grid);
  solver.setMaxIterations(1e6);
  solver.set_solver(SOR_MPI);
  solver.setRelaxationParameter(1.9);
  //solver.set_exactsol_fun(exactsol);
  //solver.exact_solution();

  solver.Solve();

  solver.WriteSolution("sol.dat");


  ierr = MPI_Finalize();
  return 0;
}
