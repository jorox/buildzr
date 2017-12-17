#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "unitcell.h"
#include "atom.h"
#include <vector>
#include "Eigen/Dense"
#include "Eigen/StdVector"

class Crystal{

 public:

  Crystal();
  Crystal(int, const Eigen::Matrix<float, 3, 2>& , UnitCell* );

  int build(std::vector<Atom>&) const;
  void print_me() const;
  int number_of_atoms() const;

  UnitCell* uc;
  Eigen::Matrix<float, 3, 2> tiles;
  int _id;
};

#endif
