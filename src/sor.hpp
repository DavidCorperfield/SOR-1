#ifndef SOR_HPP
#define SOR_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include "poissoneq.hpp"
#include "grid.hpp"
#include "boundaryconditions.hpp"

using namespace std;

class SOR{

private:

  double** solution;
  double** solution_old;
  Grid* grid;
  BoundaryConditions* bc;
  PoissonEq* pe;
  double tolerance;
  long unsigned int itmax;
  double omega;
  double InitGuess;
  double (*exactsol)(const double&, const double&);

public:

   SOR(PoissonEq* p, BoundaryConditions* b,Grid* g){
     tolerance = 1e-5;
     itmax = 1000;
     omega = 1.0;
     InitGuess = 0;
     pe = p;
     bc = b;
     grid = g;
     solution = new double* [grid->sizex()];
     solution_old = new double* [grid->sizex()];
     for(int i = 0; i < grid->sizex(); i++){
         solution[i] = new double [grid->sizey()];
         solution_old[i] = new double [grid->sizey()];
       }
     SetInitialGuess();
     setBC();
   }

   void setTolerance(double tol){
     tolerance = tol;
   }

   void setMaxIterations(int x){
     itmax = x;
   }

   void setRelaxationParameter(double x){
     omega = x;
   }

   void setBC(){
     for(int i = 0; i < grid->sizex(); i++){
         for(int j = 0; j < grid->sizey(); j++){

             if(i == 0){
                 solution_old[i][j] = bc->leftBC(grid->grid[i][j].x,grid->grid[i][j].y);
                 solution[i][j] = solution_old[i][j];
               }
             if(i == grid->sizex()-1){
                 solution_old[i][j] = bc->rightBC(grid->grid[i][j].x,grid->grid[i][j].y);
                 solution[i][j] = solution_old[i][j];
               }
             if(j == 0){
                 solution_old[i][j] = bc->bottomBC(grid->grid[i][j].x,grid->grid[i][j].y);
                 solution[i][j] = solution_old[i][j];
               }
             if(j == grid->sizey()-1){
                 solution_old[i][j] = bc->topBC(grid->grid[i][j].x,grid->grid[i][j].y);
                 solution[i][j] = solution_old[i][j];
               }

           }
       }

   }

   void SetInitialGuess(){
     for(int i = 1; i < grid->sizex()-1; i++){
         for(int j = 1; j < grid->sizey()-1; j++){
             solution_old[i][j] = 0.0;
           }
       }
   }


   void solve_jacobi_omp(){

     double a = pe->mCoeffofUxx/(2*pe->mCoeffofUxx + 2*pe->mCoeffofUyy*pow(grid->beta,2));
     double b = pe->mCoeffofUyy*pow(grid->beta,2)/(2*pe->mCoeffofUxx + 2*pe->mCoeffofUyy*pow(grid->beta,2));

     for(long unsigned int it = 0; it < itmax; it++){

         int I = grid->sizex()-1;
         int J = grid->sizey()-1;
         double sum = 0;
         int i,j;

#pragma omp parallel reduction(+:sum)
         {
#pragma omp for private(i,j)
         for(i = 1; i < I; i++){
             for(j = 1; j < J; j++){

                 double rhs = pe->mpRhsFunc(grid->grid[i][j].x,grid->grid[i][j].y)*pow(grid->dx,2);

                 /* simple jacobi scheme for laplace equation */
                 //solution[i][j] = 0.25*(solution_old[i+1][j] + solution_old[i-1][j]
                 //    + solution_old[i][j+1] + solution_old[i][j-1]);

                 /* generalied jacobi for poisson eq */
                 solution[i][j] = a*(solution_old[i+1][j] + solution_old[i-1][j])
                     + b*(solution_old[i][j+1] + solution_old[i][j-1]) - rhs;

                 /* SOR --> runs only sequentially */
                 //solution[i][j] = (1-omega)*solution_old[i][j] + omega/(2*(1+pow(grid->beta,2)))
                 //    *(solution_old[i+1][j] + solution[i-1][j] + pow(grid->beta,2)*(solution_old[i][j+1]+solution[i][j-1]));

                 sum += pow(fabs(solution[i][j] - solution_old[i][j]),2);
               }
           }
         }
         //double norm = computenorm();
         double norm = sqrt(sum);
         //cout << it << " : Residual norm = " << norm << endl;
         if(norm < tolerance){
             cout << it+1 << " : Residual norm = " << norm << endl;
             for(int p = 1; p < I; p++){
                for(int q = 1; q < J; q++){
                   solution_old[p][q] = solution[p][q];
                }
             }
             break;
           }

#pragma omp parallel
{
#pragma omp for private(i,j)
        for(i = 1; i < I; i++){
           for(j = 1; j < J; j++){
              solution_old[i][j] = solution[i][j];
           }
        }
}



       }

   }


   void solve_sor_omp(){

     double coeff_rhs = pow(grid->dx,2)/pe->mCoeffofUxx;
     double beta = pe->mCoeffofUyy/pe->mCoeffofUxx*pow(grid->beta,2);

     for(long unsigned int it = 0; it < itmax; it++){

         int I = grid->sizex()-1;
         int J = grid->sizey()-1;
         double err_sum, old_val, norm;
         int i,j;
         err_sum = 0;
#pragma omp parallel reduction(+:err_sum)
{
#pragma omp for private(i,j,old_val)
         for(i = 1; i < I; i++){
             for(j = 1; j < J; j++){

                 /* red/odd-cells */
                 if((i+j)%2 != 0){
                     old_val = solution_old[i][j];
                     double rhs = pe->mpRhsFunc(grid->grid[i][j].x,grid->grid[i][j].y);
                     solution_old[i][j] = (1-omega)*solution_old[i][j] + omega/(2*(1+beta))
                         *(solution[i+1][j] + solution[i-1][j] + beta*(solution[i][j+1]+solution[i][j-1])
                         - coeff_rhs*rhs);

                     err_sum += pow((old_val - solution_old[i][j]),2);
                   }
                 if((i+j)%2 == 0){
                     old_val = solution[i][j];
                     double rhs = pe->mpRhsFunc(grid->grid[i][j].x,grid->grid[i][j].y);
                     solution[i][j] = (1-omega)*solution[i][j] + omega/(2*(1+beta))
                         *(solution_old[i+1][j] + solution_old[i-1][j]
                         + beta*(solution_old[i][j+1]+solution_old[i][j-1]) - coeff_rhs*rhs);

                     err_sum += pow((old_val - solution[i][j]),2);
                   }
               }
           }
} // end omp parallel

         norm = sqrt(err_sum);
         //cout << it+1 <<  " : Norm = " << norm << endl;
         if(norm < tolerance){
             cout << it+1 <<  " : Norm = " << norm << endl;
             for(i = 1; i < I; i++){
                 for(j = 1; j < J; j++){
                     if((i+j)%2 != 0){
                         solution[i][j] = solution_old[i][j];
                       }
                   }
               }
             break;
           }


       }

   }



   double computenorm(){
     double sum = 0;
     for(int i = 1; i < grid->sizex()-1; i++){
         for(int j = 1; j < grid->sizey()-1; j++){

             sum += pow(fabs(solution[i][j] - solution_old[i][j]),2);

           }
       }

     return sqrt(sum);
   }

   void set_exactsol_fun(double (*pFunc)(const double&, const double&)){
     exactsol = pFunc;
   }

   void exact_solution(){
     for(int i = 1; i < grid->sizex()-1; i++){
         for(int j = 1; j < grid->sizey()-1; j++){
             solution[i][j] = exactsol(grid->grid[i][j].x,grid->grid[i][j].y);
           }
       }
   }



   void WriteSolution(string filename){
     ofstream file(filename);

     for(int i = 0; i < grid->sizex(); i++){
         for(int j = 0; j < grid->sizey(); j++){
             file << solution[i][j] << "\t";
           }
         file << endl;
       }
   }



   ~SOR(){
     for(int i = 0; i < grid->sizex(); ++i){
         delete[] solution[i];
         delete[] solution_old[i];
       }
     delete[] solution;
     delete[] solution_old;
   }

};

#endif // SOR_HPP
