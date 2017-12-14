#ifndef ENKI_H
#define ENKI_H

#include "Eigen/Dense"
#include "unitcell.h"

namespace enki{
  static constexpr double _EPS_ = 1.e-4;
  int ratio_of_atoms( const UnitCell& old, const UnitCell& shiny);
  bool compare_eps(Eigen::Vector3f, Eigen::Vector3f);
  void transform (const UnitCell& old, Eigen::Matrix3f& miller, UnitCell& shiny);

};

#endif
