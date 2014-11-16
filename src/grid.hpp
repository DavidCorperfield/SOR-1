#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>
#include "mpi.h"

using namespace std;

/* Error Check */
#define CHKERRQ(n) if(n != MPI_SUCCESS) printf("Error!! Check line number : %d\n",__LINE__)

class coordinate{
  friend class grid;
public:
  double x,y;
};


class Grid{
  friend class SOR;
  friend class BoundaryConditions;
private:

  int IMAX, JMAX;
  int lIMAX, lJMAX;
  int i_lb, i_ub, j_lb, j_ub;    // lower and upper bound in i,j dir
  double dx, dy;
  double beta;
  int ierr;
  vector<vector<coordinate>> grid;

public:

  Grid(){
    IMAX = 0; JMAX = 0;
    lIMAX = 0; lJMAX = 0;
    i_lb = 0; i_ub = 0; j_lb = 0; j_ub = 0;
    dx = 0;
    dy = 0;
    beta = 0;
  }

  Grid(const int I, const int J, const PoissonEq& pe){
    SetGridSize(I,J);
    SetEquation(pe);
  }


  void SetGridSize(const int I, const int J){
    IMAX = I;
    JMAX = J;
    int rank, size;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size); CHKERRQ(ierr);

    /* compute local array size */
    lIMAX = IMAX;
    lJMAX = JMAX/size;
    /* if rows are not exactly divisible by no of proc, add to last proc */
    if(rank == size - 1){
      lJMAX += JMAX%size;
    }

    /* compute lower and upper bound */
    i_lb = 0; i_ub = IMAX;
    j_lb = rank*(JMAX/size);
    j_ub = (rank+1)*(JMAX/size);
    if(rank == size - 1){
      j_ub += JMAX%size;
    }

  }

  void SetEquation(const PoissonEq& pe){

    assert(IMAX != 0 || JMAX != 0);
    dx = fabs(pe.Xmax-pe.Xmin)/(IMAX-1);
    dy = fabs(pe.Ymax-pe.Ymin)/(JMAX-1);
    beta = dx/dy;

    grid.resize(lIMAX,vector<coordinate>(lJMAX));
    for(int i = 0; i < lIMAX; i++){
      for(int j = 0; j < lJMAX; j++){
        coordinate coord;
        coord.x = pe.Xmin + (i+i_lb)*dx;
        coord.y = pe.Ymin + (j+j_lb)*dy;
        grid[i][j] = coord;
      }
    }

  }

  int sizex(){
    return IMAX;
  }

  int sizey(){
    return JMAX;
  }

  int lsizex(){
    return lIMAX;
  }

  int lsizey(){
    return lJMAX;
  }


  ~Grid(){
    // free grid memory
    vector<vector<coordinate>>().swap(grid);
  }

};

#endif // GRID_HPP
