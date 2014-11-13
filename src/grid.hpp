#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <cmath>
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
  int i_lb, i_ub, j_lb, j_ub;    // lower and upper bound in i,j dir
  double dx, dy;
  double beta;
  int ierr;
  vector<vector<coordinate>> grid;

public:

  Grid(){
    IMAX = 0; JMAX = 0;
    i_lb = 0, i_ub = 0, j_lb = 0, j_ub = 0;
    dx = 0;
    dy = 0;
    beta = 0;
  }

  Grid(const int I, const int J, const PoissonEq& pe){
    IMAX = I;
    JMAX = J;
    dx = fabs(pe.Xmax-pe.Xmin)/(IMAX-1);
    dy = fabs(pe.Ymax-pe.Ymin)/(JMAX-1);
    beta = dx/dy;

    grid.resize(IMAX,vector<coordinate>(JMAX));
    for(int i = 0; i < IMAX; i++){
        for(int j = 0; j < JMAX; j++){
            coordinate coord;
            coord.x = pe.Xmin + i*dx;
            coord.y = pe.Ymin + j*dy;
            grid[i][j] = coord;
          }
      }
  }


  void SetGridSize(const int I, const int J){
    IMAX = I;
    JMAX = J;
  }

  int sizex(){
    return IMAX;
  }

  int sizey(){
    return JMAX;
  }



  ~Grid(){
    // free grid memory
    vector<vector<coordinate>>().swap(grid);
  }

};

#endif // GRID_HPP
