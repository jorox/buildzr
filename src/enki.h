#ifndef ENKI_H
#define ENKI_H

#include "Eigen/Dense"
#include "unitcell.h"

namespace enki{
  int ratio_of_atoms( const UnitCell& old, const UnitCell& shiny);
  void transform (const UnitCell& old, Eigen::Matrix3f& miller, UnitCell& shiny);

};

#endif
