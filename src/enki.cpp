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

void enki::transform (const UnitCell& old, const Eigen::Matrix3f& miller, UnitCell& shiny){
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
      dummy(j) -= floor(dummy(j)+_EPS_); // return to central image
      if (fabs(dummy(j)) < 0.0001){dummy(j) = 0.0;} //set any number less that 1e-5 to zero
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


void enki::get_box_vectors (const UnitCell &cell,
                            const Eigen::Matrix<float,3,2>& tiles,
                            Eigen::Matrix<float, 3, 4> &boxMatrix){

  /* cast to float */
  //Eigen::Matrix<float,3,2> tiles = tilings.cast<float>();
  /* The origin */
  boxMatrix.block<3,1>(0,3) << cell.ucv * tiles.block<3,1>(0,0);

  /* The scale up matrix is a diagonal matrix made from the repeat integers */
  Eigen::Matrix3f scaleUp = Eigen::Matrix3f::Identity();
  for (int i = 0; i < 3; ++i){
    scaleUp(i,i) =  tiles(i,1) - tiles(i,0);
  }
  /*
   The column vectors are now the box vectors.
   Note that this may not be aligned with the x-axis depending on the unit cell
   To make it LAMMPS compatible i.e. Upper Triangular the unit cell must be rotated
  */
  boxMatrix.block<3,3>(0,0) = cell.ucv * scaleUp;

}

/**
   Rotates a unit cells' vectors so that [a1] is aligned with <x>
**/
void jallisa (UnitCell& cell){
  Eigen::PartialPivLU<Eigen::Matrix3f> luCell( cell.ucv );
  cell.ucv = luCell.matrixLU().triangularView<Eigen::Upper>();
}

/**
   wrap atoms outside the box
**/
void enki::wrap_atoms_box(std::vector<Atom>& atoms,
                          const Eigen::Matrix<float, 3, 4> &box){
  // calculate the box vectors

  Eigen::Matrix3f T = box.block<3,3>(0,0); //transformation matrix
  Eigen::Matrix3f invT = T.inverse();

  Eigen::Vector3f origin = box.block<3,1>(0,3);

  Eigen::Vector3f rAtom;
  for (int i = 0; i < atoms.size(); ++i){
    rAtom = atoms[i].coords - origin; // calculate distance from origin (llc)
    rAtom = invT * rAtom; // 0 <= rAtom(i) < 1 // transform to box space
    for (int j = 0; j < 3; ++j){
      rAtom(j) -= floor(rAtom(j));  // return to central image
    }
    rAtom = T * rAtom; // return to standard space and shift
    rAtom = rAtom + origin;
    atoms[i].coords = rAtom;
  }
}
