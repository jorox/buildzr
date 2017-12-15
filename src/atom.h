#ifndef ATOM_H
#define ATOM_H

#include "Eigen/Dense"

using namespace Eigen;

class Atom{

 public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Atom();
  Atom( int, Vector4i&, Vector3f&);

  int _id;
  Vector4i ucc;
  Vector3f coords;

};

#endif
