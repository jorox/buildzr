#include "enki.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include "Eigen/StdVector"

int enki::ratio_of_atoms( const UnitCell& old, const UnitCell& shiny){
  return floor( shiny.volume() / old.volume() );
}

void enki::transform (const UnitCell& old, Eigen::Matrix3f& miller, UnitCell& shiny){
  // calculate the miller vectors in the standard basis
  shiny.ucv << old.ucv * miller;

  // make sure volume >= 1
  int roa = ratio_of_atoms( old, shiny);
  if (roa < 1){
    std::cout << "Error enki::transform() Miller indices chosen do not represent a valid transformation" << std::endl;
  }
  else{
    printf ( "... %i atoms in new basis\n", roa * old.number_of_atoms());
  }

  // determine the tiles needed to find atoms in new basis
  Eigen::Matrix<float, 3, 2> repeats;
  //std::cout << miller.colwise().minCoeff().transpose() << std::endl;;
  repeats.block<3,1>(0,0) << miller.colwise().minCoeff().transpose();
  repeats.block<3,1>(0,1) << miller.colwise().maxCoeff().transpose();

  // generate the atoms i.e. get coordinates in standard basis (Angstroms)
  std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f>> atomsOld;
  old.generate(repeats, atomsOld);
  printf ( "    generated %i atoms using old basis\n", atomsOld.size() );

  // transform atoms from standard basis to new-basis (fractional coordinates)
  std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f>> atomsShiny;
  Eigen::Vector3f dummy;
  for ( int i = 0; i < atomsOld.size(); ++i ){
    dummy = shiny.ucv.inverse() * atomsOld[i];
    for (int j = 0; j < 3; ++j){
      if (fabs(dummy(j)) < 0.00001){dummy(j) = 0.0;} //set any number less that 1e-5 to zero
      dummy(j) -= floor(dummy(j)); // return to central image
    }
    atomsShiny.push_back(dummy); // consider changing to modify "atomsOld"
  }

  if ( atomsShiny.size() != roa * old.number_of_atoms() ){
    std::cout << "! WARNING: ENKI GENERATED WRONG NUMBER OF ATOMS!" << std::endl;
  }
  std::cout << "Enki:: New-basis Atoms:" << std::endl;
  for ( int i = 0; i < atomsShiny.size(); ++i ){
    std::cout << std::endl;
    std::cout << atomsShiny[i] << std::endl;
  }


}
