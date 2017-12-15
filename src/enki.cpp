#include "enki.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include "Eigen/StdVector"

int enki::ratio_of_atoms( const UnitCell& old, const UnitCell& shiny){
  return floor( shiny.volume() / old.volume() );
}

bool enki::mysearch(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > &atoms,
                    Eigen::Vector3f& v2){
  for (int i = 0; i < atoms.size(); ++i){
    if ( (v2-atoms[i]).norm() < 5 * _EPS_ ){ return true; }
  }
  return false;
}

void enki::transform (const UnitCell& old, Eigen::Matrix3f& miller, UnitCell& shiny){
  std::cout << "... ENKI: changing unit cell, ";

  // calculate the miller vectors in the standard basis
  shiny.ucv << old.ucv * miller;

  // make sure volume >= 1
  int roa = ratio_of_atoms( old, shiny);
  if (roa < 1){
    std::cout << "Error enki::transform() Miller indices chosen do not represent a valid transformation" << std::endl;
  }
  else{
    printf ( "expecting %i atoms in the new basis based on volume analysis\n", roa * old.number_of_atoms());
  }

  // determine the tiles needed to find atoms in new basis
  Eigen::Matrix<float, 3, 2> repeats;
  //std::cout << miller.colwise().minCoeff().transpose() << std::endl;;
  repeats.block<3,1>(0,0) << miller.rowwise().minCoeff(); //choose the minimum h, k, l
  repeats.block<3,1>(0,1) << miller.rowwise().maxCoeff(); //choose the maximum h, k, l

  // take the integer Miller indices for the repeats
  for (int i = 0; i < 3; ++i){ repeats(i,0) = floor(repeats(i,0))-1.0; }
  for (int i = 0; i < 3; ++i){ repeats(i,1) = ceil(repeats(i,1))+1.0; }

  // generate the atoms i.e. get coordinates in standard basis (Angstroms)
  std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f>> atomsOld;
  old.generate(repeats, atomsOld);
  printf ( "    generated %i atoms using old basis\n", atomsOld.size() );

  // transform atoms from standard basis to new-basis (fractional coordinates)
  std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f> > atomsShiny;
  Eigen::Vector3f dummy;
  for ( int i = 0; i < atomsOld.size(); ++i ){
    dummy = shiny.ucv.inverse() * atomsOld[i];
    for (int j = 0; j < 3; ++j){
      if (fabs(dummy(j)) < 0.00001){dummy(j) = 0.0;} //set any number less that 1e-5 to zero
      dummy(j) -= floor(dummy(j)); // return to central image
    }
    atomsShiny.push_back(dummy);
  }

  // check for duplicates
  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > uniqueAtomsShiny;
  for (int i = 0; i < atomsShiny.size(); ++i){
    if ( !enki::mysearch(uniqueAtomsShiny, atomsShiny[i]) ){
    //if (true){
      uniqueAtomsShiny.push_back(atomsShiny[i]);
      std::cout << "    pushing 1 unique i = " << i << "  " << atomsShiny[i].transpose() << std::endl;
    }
  }

  // make sure you have the right number of atoms
  if ( uniqueAtomsShiny.size() != roa * old.number_of_atoms() ){
    std::cout << "ENKI; ERROR - cannot find the expected number of atoms" << std::endl;
  }
  shiny.basis.resize(3, uniqueAtomsShiny.size());
  for (int i = 0; i < uniqueAtomsShiny.size(); ++i){
    shiny.basis.block<3,1>(0,i) << uniqueAtomsShiny[i];
  }

}
