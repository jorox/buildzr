#include "Eigen/Dense"
#include "unitcell.h"
#include <cmath>
#include "enki.h"


int main(int argc, char** argv){
  double _PI_ = 3.1415926535897;
  double cos120 = cos(120./180.*_PI_);
  double sin120 = sin(120./180.*_PI_);
  // Primitive
  Eigen::Matrix3f hcpP_base;
  Eigen::Matrix<float, 3, Eigen::Dynamic> hcpP_mot(3,2);

  hcpP_base << 1.0, cos120, 0.0,
    0.0, sin120, 0.0,
    0.0, 0.0, 1.0;

  hcpP_mot << 0.0, 0.666667,
    0.0, 0.333333,
    0.0, 0.5;

  UnitCell hcpP(hcpP_base, hcpP_mot);
  hcpP.print_me();

  // Orthogonal
  Eigen::Matrix3f miller_hcpO_base;
  Eigen::Matrix<float, 3, Eigen::Dynamic> motif_O;
  miller_hcpO_base << 1., 1., 0,
    0, 2, 0,
    0, 0, 1;

  UnitCell tmp(hcpP_base, motif_O);
  tmp.print_me();

  enki::transform(hcpP, miller_hcpO_base, tmp);
  return 0;

}
