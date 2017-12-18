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

int enki::create_perfect( const UnitCell& cell,
                    const Eigen::Matrix<float,3,2>& Nx,
                    std::vector<Atom>& atoms,
                    Eigen::Matrix<float,3,4>& box){

  UnitCell cell0(cell); // create a temporary unit cell
  Crystal xtal(0, Nx, &cell0);
  int Natoms = xtal.build(atoms);
  enki::get_box_vectors(cell, Nx, box);
  return Natoms;

}

int enki::create_edge_xz(const UnitCell& cell,
                   const Eigen::Matrix<float,3,2>& Ntiles,
                   std::vector<Atom>& atoms1,
                   Eigen::Matrix<float,3,4>& box1,
                   std::vector<Atom>& atoms2,
                   Eigen::Matrix<float,3,4>& box2){

  //
  int Nx = Ntiles(0,1) - Ntiles(0,0);
  int Ny = Ntiles(1,1) - Ntiles(1,0);
  int Nz = Ntiles(2,1) - Ntiles(2,0);
  float b = cell.ucv.block<3, 1>(0, 0).norm();

  // find the mid-plane along the Z-direction
  int midz = Ntiles(2,1) + Ntiles(2,0);
  if (midz % 2 != 0){
    midz++;
  }
  midz = (int)midz / 2;

  // Create Lambda Mu unit cells
  UnitCell ucLambda(cell);
  UnitCell ucMu(cell);
  ucLambda.ucv(0, 0) -= 0.5*b/Nx;
  ucMu.ucv(0,0) += 0.5*b/Nx;

  // Repeats
  Eigen::Matrix<float,3,2> Nlambda = Ntiles;
  Eigen::Matrix<float,3,2> Nmu = Ntiles;
  Nmu(0,1)--;
  Nmu(2,1) = midz;
  Nlambda(2,0) = midz;

  // Crystals
  Crystal mu(0, Nmu, &ucMu);
  Crystal lambda(1, Nlambda, &ucLambda);

  // Atoms and boxes
  int N1 = mu.build(atoms1);
  int N2 = lambda.build(atoms2);
  enki::get_box_vectors(ucMu, Nmu, box1);
  enki::get_box_vectors(ucLambda, Nlambda, box2);

  return N1+N2;

}

int enki::create_screw_xz( const UnitCell& cell,
                     const Eigen::Matrix<float,3,2>& Nx,
                     std::vector<Atom>& atoms1,
                     Eigen::Matrix<float,3,4>& box1,
                     std::vector<Atom>& atoms2,
                     Eigen::Matrix<float,3,4>& box2){

  float b = cell.ucv.block<3,1>(0,0).norm();
  float dy = cell.ucv.block<3,1>(0,1).norm();
  float Ny = Nx(1,1) - Nx(1,0);
  float Nz = Nx(2,1) - Nx(2,0);

  Eigen::Matrix3f ucShift;
  ucShift << 0., -0.5*b/Ny, 0.,
    0., 0., 0.,
    0., 0., 0.;

  UnitCell lambdaUC(cell);
  lambdaUC.ucv += ucShift;
  //std::cout << "ddd lambda_uc=\n" << lambdaUC.ucv << std::endl;
  //ucShift(0,1) *= -1.;
  UnitCell muUC(cell);
  muUC.ucv -= ucShift;
  //std::cout << "ddd mu_uc=\n" << muUC.ucv << std::endl;
  //std::cout << "ddd uc=\n" << cell.ucv << std::endl;

  Eigen::Matrix<float,3,2> Nlambda = Nx;
  Eigen::Matrix<float,3,2> Nmu = Nx;
  Nmu(2,1) = Nmu(2,0) + (int)Nz/2;
  Nlambda(2,0) = Nmu(2,1);

  //std::cout << "ddd Nmu=\n" << Nmu << std::endl;
  //std::cout << "ddd Nlambda=\n" << Nlambda << std::endl;

  Crystal lambda (0, Nlambda, &lambdaUC);
  Crystal mu (1, Nmu, &muUC);

  int N1 = mu.build(atoms1);
  int N2 = lambda.build(atoms2);

  enki::get_box_vectors(muUC, Nmu, box1);
  enki::get_box_vectors(lambdaUC, Nlambda, box2);

  box2(0,1) = box1(0,1);

  return N1+N2;


}

/**
   Create a SIA loop by adding atoms to the list of atoms.
   @cell is the unit cell used to create the crystal
   @tiles is the number of tilings of the unit cell
   @loopMiller 3x2 matrix representing the smallest element of the loop. Normally these will be integers
   @loopTiles a 2x2 matrix defnining how many tiles of the loop to create. These are in loop-space
   @siaV a vector representing the intersitial atom direction in loop-space
   @atoms the vector of append to
 **/

int enki::create_sia_loop ( const UnitCell& cell,
                      const Eigen::Matrix<float,3,2>& tiles,
                      const Eigen::Matrix<float,3,3>& loopMiller,
                      const Eigen::Matrix<float,3,2>& loopTiles,
                      const Eigen::Vector3f& sia,
                      std::vector<Atom>& atoms,
                      Eigen::Matrix<float,3,4>& box ){

  // First we need to find the unit cell of the loop
  UnitCell loopCell;
  // loop contains the vectors (and base atoms) of the loop in real-space
  enki::transform(cell, loopMiller, loopCell);
  float eps = 0.001;

  std::cout << "ddd new loopCell = " << std::endl;
  loopCell.print_me();

  for (int i = 0; i < loopCell.number_of_atoms(); ++i){
    loopCell.basis.block<3,1>(0,i) = loopCell.basis.block<3,1>(0,i) + (eps * sia);
  }

  std::cout << "ddd new loopCell = " << std::endl;
  loopCell.print_me();

  // ensure that at least one of the tilings has a length of 1
  /*
  bool isLoop = false;
  for (int i = 0; i < 3; ++i){
    if ( loopTiles(i,1) - loopTiles(i,0) < 1.01 ) { isLoop = true;}
  }
  if (false == isLoop) {
    std::cout << "ERROR: Cannot create loop. Must have one side length equal to 1" << std::cout;
    return 0;
  }
  */
  std::cout << "ddd loopTIles = " << std::endl << loopTiles << std::endl;
  Crystal loop(-1, loopTiles, &loopCell);
  int nLoopAtoms = loop.build(atoms);
  enki::get_box_vectors (loopCell, loopTiles, box);

  // Now loop over all atoms and find inside the loop
  /*
  int nLoopAtoms = 0;
  int oldN = atoms.size();
  Eigen::Matrix3f invLoop = loop.block<3,3>(0,0).inverse();
  Eigen::Vector3f scAtom, rAtom;

  for (int ia = 0; ia < oldN ; ++ia){
    scAtom = invLoop * atoms[ia].coords; // calculate coorindates in crystal-space
    scAtom = scAtom - loop.block<3,1>(0,3); // shift to loop origin
    std::cout << "ddd test atom = " << scAtom.transpose() << std::endl;
    for (int ix = 0; ix < 3; ++ix){ if (abs(scAtom(ix)) > 0){ continue; } }
    std::cout << "ddd loop atom = " << scAtom.transpose() << std::endl;
    scAtom = scAtom + loop.block<3,1>(0,3); // shift back to position
    scAtom = scAtom + (eps * sia); // add a small shift SIA
    rAtom = loop.block<3,3>(0,0) * scAtom; // calculate coordinates in real space
    atoms.push_back(Atom(oldN + nLoopAtoms, scAtom, rAtom)); // add SIA
    nLoopAtoms++;
  }
  */

  return nLoopAtoms;

}

/**
   Write an atom list to a LAMMPS data file
   Only orthogonal box, 1 atom type, 1 header, and 1 orientation defined for the moment
 **/
int enki::write_lammps_data_file( std::FILE* fstream,
                            std::vector<Atom>& atomList,
                            Eigen::Matrix<float,3,4>& boxLims ){

  // count number of valid atoms
  int natoms = atomList.size();
  for (int i = 0; i < atomList.size(); ++i){
    if (atomList[i]._id == -1){ natoms--; }
  }

  // make sure the box is upper triangular
  for (int i = 1; i < 3; ++i){
    for (int j = i-1; j < i; ++j){
      if (boxLims(i,j) > 1.e-4){
        std::cout << "ERROR: Crystal is not well oriented, no atoms written" << std::endl;
        return 0;
      }
    }
  }

  // write the header
  fprintf ( fstream, "#ver .1, HCP-Zr [11-20], [0001], [-1100], metal units\n");
  fprintf ( fstream, "    %i atoms\n", natoms);
  fprintf ( fstream, "    1 atom types\n\n");
  fprintf ( fstream, "%1.6f %1.6f xlo xhi\n", boxLims(0,3), boxLims(0,0)+boxLims(0,3));
  fprintf ( fstream, "%1.6f %1.6f ylo yhi\n", boxLims(1,3), boxLims(1,1)+boxLims(1,3));
  fprintf ( fstream, "%1.6f %1.6f zlo zhi\n", boxLims(2,3), boxLims(2,2)+boxLims(2,3));
  fprintf ( fstream, "%1.6f %1.6f %1.6f xy xz yz\n\n", boxLims(0,1), boxLims(0,2), boxLims(1,2));
  fprintf ( fstream, "Masses\n\n");
  fprintf ( fstream, "    1 91.224\n\n");
  fprintf ( fstream, "Atoms\n");

  int iatom = 0;
  for ( int i=0; i<atomList.size(); ++i){
    if (atomList[i]._id == -1){ continue; }
    iatom++;
    fprintf ( fstream, "\n\t%i\t%i\t %1.6f %1.6f %1.6f", iatom, 1,
              atomList[i].coords(0), atomList[i].coords(1), atomList[i].coords(2));
  }

  return iatom;

}
