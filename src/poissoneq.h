#ifndef POISSONEQ_H
#define POISSONEQ_H

#include <iostream>

using namespace std;

class PoissonEq
{
  friend class SOR;
  friend class Grid;
private:
  double mCoeffofUxx;
  double mCoeffofUyy;
  double (*mpRhsFunc)(double x, double y);
  double Xmin, Xmax, Ymin, Ymax;

public:
  PoissonEq(double coeffUxx, double coeffUyy, double (*righthandside)(double,double),
            double xmin, double xmax, double ymin, double ymax){

    mCoeffofUxx = coeffUxx;
    mCoeffofUyy = coeffUyy;
    mpRhsFunc = righthandside;
    Xmin = xmin;
    Xmax = xmax;
    Ymin = ymin;
    Ymax = ymax;
  }

};

#endif // POISSONEQ_H
