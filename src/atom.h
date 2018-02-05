#ifndef ATOM_H
#define ATOM_H

#include "Eigen/Dense"

using namespace Eigen;

class Atom{

 public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Atom();
  Atom( int,
        int,
        const Vector3d& vscaled,
        const Vector3d& vreal);

  int _id;
  Vector3d ucc;
  Vector3d coords;
  int _type;

};

#endif
