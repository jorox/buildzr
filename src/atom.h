#ifndef ATOM_H
#define ATOM_H

#include "Eigen/Dense"

using namespace Eigen;

class Atom{

 public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Atom();
  Atom( int, const Vector3f& vscaled, const Vector3f& vreal);

  int _id;
  Vector3f ucc;
  Vector3f coords;

};

#endif
