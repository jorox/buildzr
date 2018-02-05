#ifndef UNITCELL_H
#define UNITCELL_H
#include "atom.h"
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <vector>
#include <string>

class UnitCell{

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    UnitCell (const Eigen::Matrix3d& unitCellVectors,
              const Eigen::Matrix<double,3,Eigen::Dynamic>& motifAtoms,
              const std::vector<int>& motifTypes);
    UnitCell ();
  /*!
     Return the volume of the unit cell
   **/
  double volume() const;
  /*!
     Return the length of the three unit cell vectors
     \param sizes a vector to store the 3 lengths
   **/
  void size(Eigen::Vector3d& sizes) const;

  /*!
    return the number of basis atoms
   */
  int number_of_atoms() const;

  /*!
    Print information about the cell to the console
   */
  void print_me() const;

  int generate(const Eigen::Matrix<int,3,2>& repeatVecs,
                std::vector<Atom>& atoms) const;

  Eigen::Matrix3d ucv; //<-! the unit cell vectors
  Eigen::Matrix<double, 3, Eigen::Dynamic> basis; //<-! the basis atoms stored as column vectors
  std::vector<int> types;
  std::string name; //<-! the name for the unit cell
};
#endif
