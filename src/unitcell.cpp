#include "unitcell.h"
#include "atom.h"
#include "Eigen/Dense"
#include <iostream>

/*!
  \brief Construct using cell vectors and basis atoms
  \param unitCellvectors the unit cell vectos as column vectors
  \param motifAtoms the basis atoms in fractional coordinates as column vectors
  \param motifTypes a vector containing the type of each basis atom
  \todo Verify that motifAtoms has the same number of columns as the rows of motifTypes
 */
UnitCell::UnitCell(const Eigen::Matrix3d& unitCellVectors,
                   const Eigen::Matrix<double, 3, Eigen::Dynamic>& motifAtoms,
                   const std::vector<int>& motifTypes){
  ucv = unitCellVectors;
  basis = motifAtoms;
  for (int i=0; i < motifTypes.size(); ++i) {types.push_back(motifTypes[i]);}
}
/*!
  \brief create an empty unit cell. Unit cell vectos are initialized to I, no basis atoms
 */
UnitCell::UnitCell (){
  ucv = Eigen::Matrix<double,3,3>::Identity();
  basis.resize(3,0);
  name = "";
}

double UnitCell::volume() const{
  return abs(ucv.determinant());
}

void UnitCell::size(Eigen::Vector3d & lengths) const{
  for (int i = 0; i < 3; ++i){
    lengths(i) = ucv.block<3, 1>(0, i).norm();
  }
}

int UnitCell::number_of_atoms() const{
  return basis.cols();
}

void UnitCell::print_me() const {
  std::cout << "unit cell (" << name << "):" << std::endl;
  std::cout << "  basis <column vectors> (A) =" << std::endl << ucv << std::endl;
  std::cout << "  motif <column vectors> (frac)" << basis.cols() << " atoms = ";
  std::cout << std::endl << basis << std::endl;
}

/*!
  \brief Generate atoms by repeating the unit cell in 3D
  \param[in] repeatVecs a 3x2 Matrix of integers that represent the lowest and highest limits
  \param[out] atoms a vector to store the atoms
  \return the number of atoms created
*/
int UnitCell::generate(const Eigen::Matrix<int, 3, 2> &repeatVecs,
                        std::vector<Atom> &atoms) const{

  Eigen::Vector3d dummy;
  Eigen::Vector3d shift;
  int count = 0;

  // order of loop indices is important !
  for (int k = repeatVecs(2,0); k < repeatVecs(2,1); ++k){
    for (int j = repeatVecs(1,0); j < repeatVecs(1,1); ++j){
      for (int i = repeatVecs(0,0); i < repeatVecs(0,1); ++i){
        for (int ib = 0; ib < this->basis.cols(); ++ib){
          shift << (double) i, (double) j, (double) k;
          shift += this->basis.block<3,1>(0,ib); // crystal-space
          dummy = this->ucv * shift; // change to real-space
          atoms.push_back(Atom(count, types[ib], shift, dummy));
          count++;
        }
      }
    }
  }

  return count;

}
