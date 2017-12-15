#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "Eigen/Dense"
#include "unitcell.h"

class Crystal{

 public:

  Crystal();
  Crystal(int, Eigen::Matrix<int, 2, 3>& tiles, UnitCell& uc);

  double lx(int ) const;

  UnitCell unitCell;
  Eigen::Matrix<int, 2, 3> tiles;
};

#endif
