#ifndef ENKI_H
#define ENKI_H
#include "unitcell.h"
#include "atom.h"
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <vector>

namespace enki{
  static constexpr double _EPS_ = 1.e-4;

  bool mysearch (const std::vector<Atom>& ,
                 const Eigen::Vector3d &);

  void transform (const UnitCell& old,
                  const Eigen::Matrix3d& miller,
                  UnitCell& shiny);

  int ratio_of_atoms( const UnitCell& old,
                      const UnitCell& shiny);

  void transform (const UnitCell& old,
                  Eigen::Matrix3d& miller,
                  UnitCell& shiny);

  void get_box_vectors (const UnitCell& cell,
                        const Eigen::Matrix<int,3,2>& tiles,
                        Eigen::Matrix<double,3,4>& boxVectors);

  void jallisa (UnitCell& cell);

  void wrap_atoms_box (std::vector<Atom>& atoms,
                       const Eigen::Matrix<double,3,4>& box);

  /**
     \brief Creates Atom::Atom using a UnitCell::UnitCell and a repition patters
     \param cell The unit cell
     \param Nx A 3-by-2 matrix indicating the position of the lowest and higest unit cell
     \param atom A vector of Atom::Atom to add to
     \param box A 3-by-4 matrix representing a box

     The method produces atoms by creating tilings and then generating atoms with those tilings
   **/
  int create_perfect( const UnitCell& cell,
                      const Eigen::Matrix<int,3,2>& Nx,
                      std::vector<Atom>& atoms,
                      Eigen::Matrix<double,3,4>& box);
  /**
     \brief Create a crystal with an edge dislocation along the <a2> direction using Bacon's method
     \param cell The perfect UnitCell::UnitCell
     \param Ntiles A 3-by-2 matrix representing the tilings to use for the dislocated crystal
     \param atoms1 A vector to store the atoms in the lower half (mu) of the crystal
     \param box1 A 3-by-4 matrix representing the box holding the lower atoms
     \param atoms2 A vector to store tha atoms in the upper half (lambda) of the crystal
     \param box2 A 3-by-4 matrix representing the box holding the top part of the atoms

     The method creates a dislocated crystal by returning two halves having different lattice spacings
     along the <a1> direction.

     \attention The method should work for non-orthogonal unit cells, but has not been tested
     \todo Test the method with non-orthogonal unit cell
  **/
int create_edge_xz(const UnitCell& cell,
                   const Eigen::Matrix<int,3,2>& Ntiles,
                   std::vector<Atom>& atoms1,
                   Eigen::Matrix<double,3,4>& box1,
                   std::vector<Atom>& atoms2,
                   Eigen::Matrix<double,3,4>& box2);

int create_screw_xz( const UnitCell& cell,
                     const Eigen::Matrix<int,3,2>& Nx,
                     std::vector<Atom>& atoms1,
                     Eigen::Matrix<double,3,4>& box1,
                     std::vector<Atom>& atoms2,
                     Eigen::Matrix<double,3,4>& box2);

int create_sia_loop ( const UnitCell& cell,
                      const Eigen::Matrix<double,3,3>& loopMiller,
                      const Eigen::Matrix<int,3,2>& loopTiles,
                      const Eigen::Vector3d& sia,
                      std::vector<Atom>& atoms,
                      Eigen::Matrix<double,3,4>& box );

int write_lammps_data_file ( std::FILE* fstream,
                             std::vector<Atom>& atomList,
                             Eigen::Matrix<double,3,4>& boxLims );

};

#endif
