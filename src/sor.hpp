#ifndef SOR_HPP
#define SOR_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include "poissoneq.hpp"
#include "grid.hpp"
#include "boundaryconditions.hpp"
#include "mpi.h"

using namespace std;

/* Error Check */
#define CHKERRQ(n) if(n != MPI_SUCCESS) printf("Error!! Check line number : %d\n",__LINE__)

class SOR{

private:

  double **solution, **solution_old;
  double *sol_data, *sol_old_data;
  Grid* grid;
  BoundaryConditions* bc;
  PoissonEq* pe;
  double tolerance;
  long unsigned int itmax;
  double omega;
  double InitGuess;
  double (*exactsol)(const double&, const double&);
  int ierr;
  int rank, size;
  int sizex, sizey;
  int xs,xm,ys,ym;


public:

  SOR(PoissonEq* p, BoundaryConditions* b,Grid* g){
    tolerance = 1e-5;
    itmax = 1000;
    omega = 1.0;
    InitGuess = 0;
    pe = p;
    bc = b;
    grid = g;
    get_rank_size();
    allocate_sol_memory();
    SetInitialGuess();
    setBC();
  }

  void get_rank_size(){
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size); CHKERRQ(ierr);
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

  void allocate_sol_memory(){

    /* get solution vec ranges -> including ghost row */
    sizex = grid->lsizex();
    if(size == 1){
      sizey = grid->lsizey();
    }else if(size > 1){
      if(rank == 0 || rank == size - 1){
        sizey = grid->lsizey()+1;
      }else{
        sizey = grid->lsizey()+2;
      }
    }
    //cout << rank << " : sizey = " << sizey << endl;

    /* allocate contiguous memory for solution vec -> including ghost row */
    solution = new double* [sizex];
    solution_old = new double* [sizex];
    sol_data = new double [sizex*sizey];
    sol_old_data = new double [sizex*sizey];
    for(int i = 0; i < sizex; i++){
      solution[i] = &sol_data[i*sizey];
      solution_old[i] = &sol_old_data[i*sizey];
    }

    /* assign starting and max index of sol vec -> ignoring ghost row */
    xs = 0;
    xm = grid->lsizex();
    if(size == 1){
      ys = 0;
      ym = grid->lsizey();
    }else if(size > 1){
      if(rank == 0){
        ys = 0;
        ym = sizey-1;
      }else if(rank == size-1){
        ys = 1;
        ym = sizey;
      }else if(rank != 0 && rank != size-1){
        ys = 1;
        ym = sizey - 1;
      }
    }
    //cout << rank << " : ys = " << ys << " \t ym = " << ym << endl;

  }

  void setBC(){
    int _i, _j, j_;
    for(int i = xs; i < xm; i++){
      for(int j = ys; j < ym; j++){

        /* index offset to map to global array */
        _i = i + grid->i_lb;
        if(size == 1){
          _j = j + grid->j_lb;
        }else if(size > 1){
          if(rank == 0){
            _j = j + grid->j_lb;
          }else if(rank != 0){
            _j = (j-1) + grid->j_lb;
          }
        }

        /* index offset to map to vec with no ghost node */
        if(rank == 0){
          j_ = j;
        }else if(rank != 0){
          j_ = j-1;
        }

        if(_i == 0){
          solution_old[i][j] = bc->leftBC(grid->grid[i][j_].x,grid->grid[i][j_].y);
          solution[i][j] = solution_old[i][j];
        }
        if(_i == grid->sizex()-1){
          solution_old[i][j] = bc->rightBC(grid->grid[i][j_].x,grid->grid[i][j_].y);
          solution[i][j] = solution_old[i][j];
        }
        if(_j == 0){
          solution_old[i][j] = bc->bottomBC(grid->grid[i][j_].x,grid->grid[i][j_].y);
          solution[i][j] = solution_old[i][j];
        }
        if(_j == grid->sizey()-1){
          solution_old[i][j] = bc->topBC(grid->grid[i][j_].x,grid->grid[i][j_].y);
          solution[i][j] = solution_old[i][j];
        }

      }
    }

  }

  void SetInitialGuess(){
    int p = 0;
    for(int i = 0; i < sizex; i++){
      for(int j = 0; j < sizey; j++){
          solution_old[i][j] = 0;
          solution[i][j] = 0;
          p++;
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



  void solve_sor_mpi(){

    double norm;
    double coeff_rhs = pow(grid->dx,2)/pe->mCoeffofUxx;
    double beta = pe->mCoeffofUyy/pe->mCoeffofUxx*pow(grid->beta,2);

    /* define datatype to send row of 2d array */
    MPI_Datatype send_row;
    MPI_Type_vector(sizex,1,sizey,MPI_DOUBLE,&send_row);
    MPI_Type_commit(&send_row);

    for(long unsigned int it = 0; it < itmax; it++){

      double err_sum = 0;

      if(size > 1){
        /* send up unless I'm at the top */
        if(rank < size-1){
          ierr = MPI_Send(&solution_old[0][sizey-2],1,send_row,rank+1,0,MPI_COMM_WORLD); CHKERRQ(ierr);
        }
        /* receive from below unless I'm at the bottom */
        if(rank > 0){
          ierr = MPI_Recv(&solution_old[0][0],1,send_row,rank-1,0,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); CHKERRQ(ierr);
        }
        /* send down upless I'm at bottom */
        if(rank > 0){
          ierr = MPI_Send(&solution_old[0][1],1,send_row,rank-1,1,MPI_COMM_WORLD); CHKERRQ(ierr);
        }
        /* receive from above unless i'm at the top */
        if(rank < size-1){
          ierr = MPI_Recv(&solution_old[0][sizey-1],1,send_row,rank+1,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); CHKERRQ(ierr);
        }
      }

      for(int i = xs; i < xm; i++){
        for(int j = ys; j < ym; j++){

          int _i,_j,j_;

          /* index offset to map to global array */
          _i = i + grid->i_lb;
          if(size == 1){
            _j = j + grid->j_lb;
          }else if(size > 1){
            if(rank == 0){
              _j = j + grid->j_lb;
            }else if(rank != 0){
              _j = (j-1) + grid->j_lb;
            }
          }
          /* index offset to map to vec with no ghost node */
          if(rank == 0){
            j_ = j;
          }else if(rank != 0){
            j_ = j-1;
          }
          double rhs = pe->mpRhsFunc(grid->grid[i][j_].x,grid->grid[i][j_].y);

          /* exclude boundary nodes where BC is present */
          if(_i > 0 && _i < grid->sizex()-1 && _j > 0 && _j < grid->sizey()-1){

            solution[i][j] = (1-omega)*solution_old[i][j] + omega/(2*(1+beta))
                *(solution_old[i+1][j] + solution[i-1][j]
                + beta*(solution_old[i][j+1]+solution[i][j-1]) - coeff_rhs*rhs);

            err_sum += pow((solution[i][j] - solution_old[i][j]),2);
            solution_old[i][j] = solution[i][j];
          }

        }
      }

      if(it%50 == 0){
        MPI_Allreduce(&err_sum,&norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        norm = sqrt(norm);
        if(rank == 0){
          cout << it << " : Residual norm = " << norm << endl;
        }
      }

      if(norm < tolerance){
        break;
      }

    } /* end jacobi iteration loop */


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
    int j_;
    for(int i = xs; i < xm; i++){
      for(int j = ys; j < ym; j++){
        /* index offset to map to vec with no ghost node */
        if(rank == 0){
          j_ = j;
        }else if(rank != 0){
          j_ = j-1;
        }
        solution[i][j] = exactsol(grid->grid[i][j_].x,grid->grid[i][j_].y);
      }
    }
  }



  void WriteSolution(string filename){

    string name = to_string(rank);
    name += filename;
    ofstream file(name);
    for(int j = ys; j < ym; j++){
      for(int i = xs; i < xm; i++){
        file << scientific << solution[i][j] << "\t";
      }
      file << endl;
    }
  }



  ~SOR(){
    delete[] solution;
    delete[] solution_old;
    delete [] sol_data;
    delete [] sol_old_data;
  }

};

#endif // SOR_HPP
