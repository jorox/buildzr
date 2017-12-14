#include "unitcell.h"
#include "Eigen/Dense"
#include <iostream>

UnitCell::UnitCell(Eigen::Matrix3f& unitCellVectors, Eigen::Matrix<float, 3, Eigen::Dynamic>& motifAtoms){
  ucv = unitCellVectors;
  basis = motifAtoms;
}

UnitCell::UnitCell (){
  ucv = Eigen::Matrix<float,3,3>::Identity();
  basis.resize(3,0);
}

double UnitCell::volume() const{
  return ucv.determinant();
}

void UnitCell::size(Eigen::Vector3f & lengths) const{
  for (int i = 0; i < 3; ++i){
    lengths(i) = ucv.block<3, 1>(0, i).norm();
  }
}

int UnitCell::number_of_atoms() const{
  return basis.cols();
}

void UnitCell::print_me() const {
  std::cout << "basis <column vectors> (A) =" << std::endl << ucv << std::endl;
  std::cout << "motif <column vectors> (frac) =" << std::endl << basis << std::endl;
}

void UnitCell::generate(const Eigen::Matrix<float,3,2> &repeatVecs, std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > &atoms) const{

  Eigen::Vector3f dummy;
  Eigen::Vector3f shift;

  // order of loop indices is important !
  for (int k = repeatVecs(2,0); k < repeatVecs(2,1); ++k){
    for (int j = repeatVecs(1,0); j < repeatVecs(1,1); ++j){
      for (int i = repeatVecs(0,0); i < repeatVecs(0,1); ++i){
        for (int ib = 0; ib < this->basis.cols(); ++ib){
          shift << i, j, k;
          dummy = this->ucv * (this->basis.block<3,1>(0,ib) + shift);
          atoms.push_back(dummy);
        }
      }
    }
  }


}
