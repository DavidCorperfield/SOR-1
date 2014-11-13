#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "poissoneq.h"
#include "boundaryconditions.h"

using namespace std;

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
  double dx, dy;
  double beta;

  vector<vector<coordinate>> grid;

public:

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

  int sizex(){
    return IMAX;
  }

  int sizey(){
    return JMAX;
  }

  double DX(){
    return dx;
  }

  double DY(){
    return dy;
  }

  double Beta(){
    return beta;
  }


  ~Grid(){
    // free grid memory
    vector<vector<coordinate>>().swap(grid);
  }

};

#endif // GRID_H
