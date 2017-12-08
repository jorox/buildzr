#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"

struct UnitCell {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Eigen::Vector3f X, Y, Z;
  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> motif;
};

struct Atom {
  int id, ix, iy, iz, ib;
  double x, y, z;
};


void id_to_uc_coordinates(int atom_id, Eigen::Vector3i& atom_uc, Eigen::Matrix<3,2>& boxuc, int nbasis)

int main(int argc, char** argv){

  if (argc != 7){ //ensure that 6 number are given -- they will be cast as integers
    std::cout << "ERROR: must supply exactly 6 integers to determine box limits" << std::endl;
    return 1;
  }

  /*
    lattice constants
   */
  double _ALAT = 3.23; // Angstrom
  double _CLAT = 5.165;

  /*
    Zirconium unit cell
    X = 1/3[11-20]
    Y = [0001]
    Z = [-1100]
    We have 4 basis atoms in this unit cell i.e. it is not the primitive; however, it is orthogonal
   */
  struct UnitCell zrstd;

  zrstd.X << _ALAT, 0.0, 0.0; // 1/3[11-20]
  zrstd.Y << 0.0, _CLAT, 0.0; // [0001]
  zrstd.Z << 0.0, 0.0, sqrt(3.0) * _ALAT; // [1-100]

  zrstd.motif.push_back(Eigen::Vector3f (0.0, 0.0, 0.0));   //A-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.0, 0.5, 1./3.)); //B-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.5, 0.0, 0.5));   //A-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.5, 0.5, 2./3.)); //B-plane

  /*
    The box is defined by the repeat vectors Nx[:,i]
    The bottom-left and top-right corners are the multiple of the unit cell vectors with the number of unit cells
    along that direction. Note that in this way the bottom-left corner might NOT be the left-most point if the
    unit cell is not orthogonal for example.
   */
  Eigen::Matrix<int, 3, 2> Nx; // box dimensions in unit cell
  Eigen::Matrix<float, 3, 3> box_lim; // box limits i.e bottom-left and top-right corners
  int Natoms = 0;

  for (int i = 1; i<7; ++i){
    Nx((i - 1)/2, (i-1)%2) = atoi(argv[i]);
  }
  // calculate the blc and trc
  box_lim.block<3,1>(0,0) = Nx(0,0) * zrstd.X + Nx(1,0) * zrstd.Y + Nx(2,0) * zrstd.Z; //lower-left corner
  box_lim.block<3,1>(0,1) = Nx(0,1) * zrstd.X + Nx(1,1) * zrstd.Y + Nx(2,1) * zrstd.Z; //upper-right corner

  // calculate the number of atoms
  Natoms = (Nx(0,1) - Nx(0,0)) * (Nx(1,1) - Nx(1,0)) * (Nx(2,1) - Nx(2,0)) * zrstd.motif.size();

  // print information on the box
  printf ( "... box dimensions in lattice units: \n     X = [%i %i]\n     Y = [%i %i]\n     Z = [%i %i]\n",
           Nx(0,0), Nx(0,1), Nx(1,0), Nx(1,1), Nx(2,0), Nx(2,1));
  printf ( "    => %i atoms\n\n", Natoms);

  printf ( "... spacings along X, Y, Z: %1.5f, %1.5f, %1.5f\n", zrstd.X.norm(), zrstd.Y.norm(),zrstd.Z.norm());

  printf ( "... box dimensions in A: \n     X = %1.4f --> %1.4f\n     Y = %1.4f --> %1.4f\n     Z = %1.4f --> %1.4f\n\n",
           box_lim(0,0), box_lim(0,1), box_lim(1,0), box_lim(1,1), box_lim(2,0), box_lim(2,1));


  printf("... building crystal");
  std::vector<Atom> atoms;

  for (int i = 0; i < Natoms; ++i){
    atoms.push_back(Atom{i, })

  }
  return 0;

}
