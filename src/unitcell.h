#ifndef UNITCELL_H
#define UNITCELL_H
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <vector>
#include "atom.h"

class UnitCell{

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    UnitCell (Eigen::Matrix3f& unitCellVectors,
              Eigen::Matrix<float, 3, Eigen::Dynamic>& motifAtoms);
    UnitCell ();
  /**
     Return the volume of the unit cell
   **/
  double volume() const;
  /**
     Return the length of the three unit cell vectors
   **/
  void size(Eigen::Vector3f& ) const;

  int number_of_atoms() const;

  void print_me() const;

  int generate(const Eigen::Matrix<float,3,2>& repeatVecs,
                std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> >& atoms) const;

  int generate(const Eigen::Matrix<float,3,2>& repeatVecs,
                std::vector<Atom>& atoms) const;

  Eigen::Matrix3f ucv;
  Eigen::Matrix<float, 3, Eigen::Dynamic> basis;
};
#endif
