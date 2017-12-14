#include "Eigen/Dense"
#include "unitcell.h"
#include <cmath>
#include "enki.h"
#include <iostream>


int main(int argc, char** argv){
  double _PI_ = 3.1415926535897;
  double cos120 = cos(120./180.*_PI_);
  double sin120 = sin(120./180.*_PI_);

  std::cout << "Transform HCP-primitive to HCP-orthogonal" << std::endl;
  // HCP-P
  Eigen::Matrix3f hcp_P_base;
  Eigen::Matrix<float, 3, Eigen::Dynamic> hcp_P_mot(3,2);
  hcp_P_base << 1.0, cos120, 0.0,
    0.0, sin120, 0.0,
    0.0, 0.0, 1.0;
  hcp_P_mot << 0.0, 0.666667,
    0.0, 0.333333,
    0.0, 0.5;
  UnitCell hcp_P(hcp_P_base, hcp_P_mot);
  std::cout << "======================================" << std::endl;
  std::cout << "===        HCP - Primitive         ===" << std::endl;
  std::cout << "======================================" << std::endl;
  hcp_P.print_me();

  // HCP-O
  Eigen::Matrix3f hcp_O_P_miller;
  hcp_O_P_miller << 1., 1., 0,
    0, 2, 0,
    0, 0, 1;

  UnitCell hcp_O;
  enki::transform(hcp_P, hcp_O_P_miller, hcp_O);
  std::cout << "======================================" << std::endl;
  std::cout << "===        HCP - Orthogonal        ===" << std::endl;
  std::cout << "======================================" << std::endl;
  hcp_O.print_me();

  // BCC cubic conventional
  Eigen::Matrix3f bcc_c_base;
  Eigen::Matrix<float,3,Eigen::Dynamic> bcc_c_motf(3,2);
  bcc_c_base << Eigen::Matrix3f::Identity();
  bcc_c_motf << 0.0, 0.5, 0.0, 0.5, 0.0, 0.5;
  UnitCell bcc_c(bcc_c_base, bcc_c_motf);

  std::cout << "======================================" << std::endl;
  std::cout << "===       BCC - Conventional       ===" << std::endl;
  std::cout << "======================================" << std::endl;
  bcc_c.print_me();

  //BCC - Disloc
  Eigen::Matrix3f bcc_D_C_miller;
  bcc_D_C_miller << .5, -1., 1.,
    .5, -1., -1.,
    .5, 2., 0.;

  UnitCell bcc_d;

  enki::transform( bcc_c, bcc_D_C_miller, bcc_d );


  std::cout << "======================================" << std::endl;
  std::cout << "===           BCC - Disloc         ===" << std::endl;
  std::cout << "======================================" << std::endl;
  bcc_d.print_me();

  std::cout << "\n... Done with Enki tests, Ba'al be praised!\n" << std::endl;
  return 0;

}
