#ifndef ENKI_H
#define ENKI_H

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "unitcell.h"
#include <vector>

namespace enki{
  static constexpr double _EPS_ = 1.e-4;
  int ratio_of_atoms( const UnitCell& old, const UnitCell& shiny);
  bool mysearch(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > & , Eigen::Vector3f &);
  void transform (const UnitCell& old, Eigen::Matrix3f& miller, UnitCell& shiny);

};

#endif
