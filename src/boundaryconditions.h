#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <iostream>
#include <vector>

class BoundaryConditions{
  friend class SOR;
private:

  bool isLeftBCset;
  bool isRightBCset;
  bool isTopBCset;
  bool isBottomBCset;

  double (*leftBC)(const double&, const double&);
  double (*rightBC)(const double&, const double&);
  double (*topBC)(const double&, const double&);
  double (*bottomBC)(const double&, const double&);

public:

  BoundaryConditions(){
    isLeftBCset = false;
    isRightBCset = false;
    isTopBCset = false;
    isBottomBCset = false;
  }

  void SetLeftBC(double (*pFunc)(const double&, const double&)){
    leftBC = pFunc;
    isLeftBCset = true;
  }

  void SetRightBC(double (*pFunc)(const double&, const double&)){
    rightBC = pFunc;
    isRightBCset = true;
  }

  void SetTopBC(double (*pFunc)(const double&, const double&)){
    topBC = pFunc;
    isTopBCset = true;
  }

  void SetBottomBC(double (*pFunc)(const double&, const double&)){
    bottomBC = pFunc;
    isBottomBCset = true;
  }


};

#endif // BOUNDARYCONDITIONS_H
