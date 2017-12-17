#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "atom.h"
#include "unitcell.h"
#include "crystal.h"
#include "enki.h"

static constexpr double _PI_ = 3.1415926535897;

/**
   Write an atom list to a LAMMPS data file
   Only orthogonal box, 1 atom type, 1 header, and 1 orientation defined for the moment
 **/
int write_lammps_data_file( std::FILE* fstream,
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

int create_perfect( const UnitCell& cell,
                    const Eigen::Matrix<float,3,2>& Nx,
                    std::vector<Atom>& atoms,
                    Eigen::Matrix<float,3,4>& box){

  UnitCell cell0(cell); // create a temporary unit cell
  Crystal xtal(0, Nx, &cell0);
  int Natoms = xtal.build(atoms);
  enki::get_box_vectors(cell, Nx, box);
  return Natoms;

}

/**
   Create an edge dislocation with Burgers vector b = x, and n = z
   Atoms are not deleted to preserve fast id-ucc mapping.
   Instead their id is changed to -1 in which case they are not written out.
 **/
int create_edge_xz(const UnitCell& cell,
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

int create_screw_xz( const UnitCell& cell,
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

int create_sia_loop ( const UnitCell& cell, const Eigen::Matrix<float,3,2>& tiles,
                      const Eigen::Matrix<float,3,2>& loop, const Eigen::Matrix<float,3,2>& loopTiles,
                      std::vector<Atom>& atoms){


}

int main(int argc, char** argv){

  if (argc < 7){ //ensure that 6 number are given -- they will be cast as integers
    std::cout << "ERROR: must supply at least 6 integers to determine box limits" << std::endl;
    return 1;
  }

  double _ALAT = 3.23; // Angstrom
  double _CLAT = 5.165;
  Eigen::Matrix<float,3,4> boxMat;
  std::vector<Atom> atoms;

  //  Zirconium unit cell
  //  X = 1/3[11-20]
  //  Y = [0001]
  //  Z = [-1100]
  //  We have 4 basis atoms in this unit cell i.e. it is not the primitive; however, it is orthogonal

  Eigen::Matrix3f zrstd_basis;
  Eigen::Matrix<float,3,Eigen::Dynamic> zrstd_motif(3,4);
  Eigen::Matrix<float, 3, 2> Nx; // box dimensions in unit cell

  /*
  // Z - basal
  zrstd_basis << _ALAT, 0.0, 0.0, // 1/3[11-20]
    0.0, sqrt(3.0) * _ALAT, 0.0, // [1-100]
    0.0, 0.0, _CLAT; // [0001]

  zrstd_motif.block<3,1>(0,0) = Eigen::Vector3f (0.0, 0.0, 0.0);   //A-plane
  zrstd_motif.block<3,1>(0,1) = Eigen::Vector3f (0.5, 0.16666667, 0.5); //B-plane
  zrstd_motif.block<3,1>(0,2) = Eigen::Vector3f (0.5, 0.5, 0.0);   //A-plane
  zrstd_motif.block<3,1>(0,3) = Eigen::Vector3f (0.0, 0.66666667, 0.5); //B-plane
  */

  // Z - prism
  zrstd_basis << _ALAT, 0.0, 0.0, // 1/3[11-20]
    0.0, _CLAT, 0.0, // [1-100]
    0.0, 0.0, sqrt(3.0) * _ALAT; // [0001]

  zrstd_motif.block<3,1>(0,0) = Eigen::Vector3f (0.0, 0.0, 0.0);   //A-plane
  zrstd_motif.block<3,1>(0,1) = Eigen::Vector3f (0.5, 0.5, 0.16666667); //B-plane
  zrstd_motif.block<3,1>(0,2) = Eigen::Vector3f (0.5, 0.0, 0.5);   //A-plane
  zrstd_motif.block<3,1>(0,3) = Eigen::Vector3f (0.0, 0.5, 0.66666667); //B-plane
  UnitCell zrstd_uc(zrstd_basis, zrstd_motif);


  for (int i = 0; i < 6; ++i){ Nx((int)i/2, i%2) = std::atof (argv[i+1]); }


  if (argc == 7){
    // print information on the box
    printf ( "... using repeats: \n     X = [%1.0f %1.0f]\n     Y = [%1.0f %1.0f]\n     Z = [%1.0f %1.0f]\n",
           Nx(0,0), Nx(0,1), Nx(1,0), Nx(1,1), Nx(2,0), Nx(2,1));

    Eigen::Vector3f box_spc;
    zrstd_uc.size(box_spc);
    printf ( "... spacings along X, Y, Z: %1.5f, %1.5f, %1.5f\n",
           box_spc(0),
           box_spc(1),
           box_spc(2));

    printf ( "... building perfect crystal\n");
    int Natoms = create_perfect(zrstd_uc, Nx, atoms, boxMat);
    printf( "    crystal built successfully\n");
    printf ( "+++ Box (Angstrom) [h1|h2|h2|origin]:\n");
    std::cout << boxMat << std::endl;
    printf ( "+++  %i atoms", Natoms);

    std::FILE * fp;
    fp = fopen ( "zr.data", "w");
    int tmp = write_lammps_data_file ( fp, atoms, boxMat );
    printf ( "... done writing %i atoms to zr.data\n", tmp );
  }


  if (argc > 7){
    Eigen::Matrix<float,3,4> boxMat1;
    std::vector<Atom> atoms1;
    int tmp;
    if (strcmp(argv[7],"e") == 0){
      tmp = create_edge_xz (zrstd_uc, Nx, atoms, boxMat, atoms1, boxMat1 );
      printf ( "+++ created edge dislocation with b=x, n=z, %i atoms \n", tmp);
    }
    else{
      tmp = create_screw_xz(zrstd_uc, Nx, atoms, boxMat, atoms1, boxMat1);
      std::cout << "ddd boxMat=\n" << boxMat << std::endl;
      std::cout << "ddd boxMat1=\n" << boxMat1 << std::endl;
      printf ( "+++ created screw dislocation with b=x, n=z,\n");
      printf ( "    %i atoms:\n", tmp);
    }

    std::FILE * fp;
    std::FILE * fp1;
    std::FILE * fpp;
    fp = fopen ( "zrMu.data", "w");
    fp1 = fopen( "zrLm.data", "w");
    fpp = fopen( "zr.data", "w");
    tmp = write_lammps_data_file ( fp, atoms, boxMat );
    printf ( "... done writing %i atoms to zrMu.data\n", tmp );
    tmp = write_lammps_data_file (fp1, atoms1, boxMat1 );
    printf ( "    done writing %i atoms to zrLm.data\n", tmp);

    std::vector<Atom> atomsFull;
    atomsFull.reserve(atoms.size()+atoms1.size());
    atomsFull.insert(atomsFull.end(), atoms.begin(), atoms.end());
    atomsFull.insert(atomsFull.end(), atoms1.begin(), atoms1.end());

    Eigen::Matrix<float, 3, 4> boxFull;
    boxFull.block<3,2>(0,0) = 0.5 * (boxMat.block<3,2>(0,0) + boxMat1.block<3,2>(0,0));
    boxFull.block<3,1>(0,2) = boxMat.block<3,1>(0,2) + boxMat1.block<3,1>(0,2);
    boxFull.block<3,1>(0,3) = boxMat.block<3,1>(0,3);
    enki::wrap_atoms_box(atomsFull, boxFull);
    tmp = write_lammps_data_file (fpp, atomsFull, boxFull );
    printf ( "    done writing %i atoms to zr.data\n", tmp);
    }

  return 0;

}
