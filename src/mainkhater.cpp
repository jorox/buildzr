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
typedef Eigen::Matrix<float,3,4> BoxMatrix;
typedef std::vector<Atom> AtomVector;

void join (const AtomVector& atomsM,
           const BoxMatrix& boxM,
           const AtomVector& atomsL,
           const BoxMatrix& boxL,
           AtomVector& atomsFull,
           BoxMatrix& boxFull){

  // first join the atom lists
  atomsFull.reserve(atomsM.size()+atomsL.size());
  atomsFull.insert(atomsFull.end(), atomsM.begin(), atomsM.end());
  atomsFull.insert(atomsFull.end(), atomsL.begin(), atomsL.end());

  // join the boxes along the Z
  boxFull.block<3,2>(0,0) = 0.5 * (boxM.block<3,2>(0,0) + boxL.block<3,2>(0,0)); // The X and Y
  boxFull.block<3,1>(0,2) = boxM.block<3,1>(0,2) + boxL.block<3,1>(0,2); // the Z-axis
  boxFull.block<3,1>(0,3) = boxM.block<3,1>(0,3); // the origin = lower box
}



int main(int argc, char** argv){

  if (argc < 9){ //ensure that 6 number are given -- they will be cast as integers
    std::cout << "ERROR: must supply at least 6 integers and 2 chars to determine box limits" << std::endl;
    return 1;
  }

  double _ALAT = 3.23; // Angstrom
  double _CLAT = 5.165;

  //  Zirconium unit cell
  //  X = 1/3[11-20]
  //  Y = [0001]
  //  Z = [-1100]
  //  We have 4 basis atoms in this unit cell i.e. it is not the primitive; however, it is orthogonal

  Eigen::Matrix3f zrstd_b_vec;
  Eigen::Matrix<float,3,Eigen::Dynamic> zrstd_b_motif(3,4);
  Eigen::Matrix3f zrstd_p_vec;
  Eigen::Matrix<float,3,Eigen::Dynamic> zrstd_p_motif(3,4);

  // Z - basal
  zrstd_b_vec << _ALAT, 0.0, 0.0, // 1/3[11-20]
    0.0, sqrt(3.0) * _ALAT, 0.0, // [1-100]
    0.0, 0.0, _CLAT; // [0001]

  zrstd_b_motif.block<3,1>(0,0) << 0.0, 0.0, 0.0;   //A-plane
  zrstd_b_motif.block<3,1>(0,1) << 0.5, 0.16666667, 0.5; //B-plane
  zrstd_b_motif.block<3,1>(0,2) << 0.5, 0.5, 0.0;   //A-plane
  zrstd_b_motif.block<3,1>(0,3) << 0.0, 0.66666667, 0.5; //B-plane

  // Z - prism
  zrstd_p_vec << _ALAT, 0.0, 0.0, // 1/3[11-20]
    0.0, _CLAT, 0.0, // [1-100]
    0.0, 0.0, sqrt(3.0) * _ALAT; // [0001]

  zrstd_p_motif.block<3,1>(0,0) << 0.0, 0.0, 0.0;   //A-plane
  zrstd_p_motif.block<3,1>(0,1) << 0.5, 0.5, 0.16666667; //B-plane
  zrstd_p_motif.block<3,1>(0,2) << 0.5, 0.0, 0.5;   //A-plane
  zrstd_p_motif.block<3,1>(0,3) << 0.0, 0.5, 0.66666667; //B-plane

  /*
    BCC dislocation unit cell
    X = 0.5[111]
    Y = [-1-12]
    Z = [1-10]
    We have 6 basis atoms in this unit cell i.e. it is not the primitive; however, it is orthogonal
  */
  double _ALAT_BCC = 2.855324;
  Eigen::Matrix3f bcc111_vec;
  Eigen::Matrix<float,3,Eigen::Dynamic> bcc111_motif(3,6);

  bcc111_vec.block<3,1>(0,0) << sqrt(3.0)/2.0, 0., 0.; // 1/2[111]
  bcc111_vec.block<3,1>(0,1) << 0.0, sqrt(6.0), 0.; // [-1-12]
  bcc111_vec.block<3,1>(0,2) << 0.0, 0.0, sqrt(2.0); // [1-10]

  bcc111_vec.block<3,1>(0,0) *= _ALAT_BCC;
  bcc111_vec.block<3,1>(0,1) *= _ALAT_BCC;
  bcc111_vec.block<3,1>(0,2) *= _ALAT_BCC;

  bcc111_motif.block<3,1>(0,0) << 0.6666667, 0.333333, 0.0;
  bcc111_motif.block<3,1>(0,1) << 0.3333333, 0.166667, 0.5;
  bcc111_motif.block<3,1>(0,2) << 0., 0., 0.;
  bcc111_motif.block<3,1>(0,3) << 0.666667, 0.8333333, 0.5;
  bcc111_motif.block<3,1>(0,4) << 0.333333, 0.6666667, 0.0;
  bcc111_motif.block<3,1>(0,5) <<  0.0, 0.5, 0.5;


  // READ COMMAND-LINE ARGUMENTS
  Eigen::Matrix<float, 3, 2> Nx; // box dimensions in unit cell
  Eigen::Matrix<float, 3, 2> Nloop; // box dimensions in unit cell
  enum mat_type {Zrb, Zrp, Fe} mat; //material: Zr or Fe
  enum crs_type {P, E, S} crs;
  enum sia_type {Y, N} sia;

  for (int i = 0; i < 6; ++i){
    Nx((int)i/2, i%2) = std::atof (argv[i+1]); // crystal parameters
  }

  printf ( "... using repeats: \n     X = [%1.0f %1.0f]\n     Y = [%1.0f %1.0f]\n     Z = [%1.0f %1.0f]\n",
           Nx(0,0), Nx(0,1), Nx(1,0), Nx(1,1), Nx(2,0), Nx(2,1));

  if (strcmp(argv[7],"Zrb")==0){mat = Zrb;}
  if (strcmp(argv[7],"Zrp")==0){mat = Zrp;}
  if (strcmp(argv[7],"Fe")==0){mat = Fe;}
  std::cout << "mat = " << mat << std::endl;
  if (strcmp(argv[8],"P")==0){crs = P;}
  if (strcmp(argv[8],"E")==0){crs = E;}
  if (strcmp(argv[8],"S")==0){crs = S;}
  std::cout << "crs = " << crs << std::endl;

  UnitCell cell;
  switch (mat){
  case Fe:{
    cell.ucv = bcc111_vec;
    cell.basis = bcc111_motif;
    break;
  }
  case Zrb:{
    cell.ucv = zrstd_b_vec;
    cell.basis = zrstd_b_motif;
    break;
  }
  case Zrp:{
    cell.ucv = zrstd_p_vec;
    cell.basis = zrstd_p_motif;
    break;
  }
  }

  Eigen::Vector3f box_spc;
  cell.size(box_spc);
  printf ( "... spacings (x,y,z): %1.5f, %1.5f, %1.5f\n",
           box_spc(0),
           box_spc(1),
           box_spc(2) );

  int Natoms;
  AtomVector atoms;
  BoxMatrix box;
  int tmp;
  BoxMatrix boxMatMu, boxMatLambda;
  AtomVector atomsMu, atomsLambda;

  switch (crs){
  case P:{
    printf ( "... building perfect crystal\n");
    Natoms = enki::create_perfect(cell, Nx, atoms, box);
    printf( "+++ crystal built successfully, %i atoms\n", Natoms);
    printf ( "    Box (Angstrom) [h1|h2|h2|origin]:\n");
    std::cout << box << std::endl;
    break;
  }

  case E:{
    printf ("... building edge-dislocated crystal\n");
    Natoms = enki::create_edge_xz (cell, Nx, atomsMu, boxMatMu, atomsLambda, boxMatLambda );
    printf ( "+++ created edge dislocation with b=x, n=z, %i atoms \n", Natoms);
    printf ( "    Box-Mu (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsMu.size());
    std::cout << "    boxMat=\n" << boxMatMu << std::endl;
    printf ( "    Box-Lambda (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsLambda.size());
    std::cout << "    boxLambda=\n" << boxMatLambda << std::endl;
    join ( atomsMu, boxMatMu, atomsLambda, boxMatLambda, atoms, box );
    break;
  }

  case S:{
    printf ("... building scrw-dislocated crystal\n");
    Natoms = enki::create_screw_xz (cell, Nx, atomsMu, boxMatMu, atomsLambda, boxMatLambda );
    printf ( "+++ created screw dislocation with b=x, n=z, %i atoms \n", Natoms);
    printf ( "    Box-Mu (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsMu.size());
    std::cout << "    boxMat=\n" << boxMatMu << std::endl;
    printf ( "    Box-Lambda (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsLambda.size());
    std::cout << "    boxLambda=\n" << boxMatLambda << std::endl;
    join ( atomsMu, boxMatMu, atomsLambda, boxMatLambda, atoms, box );
    break;
  }
  }

  if (argc > 9){ sia = Y;}
  else{ sia = N;}

  if (sia==Y){
    for (int i = 0; i < 6; ++i){
      Nloop((int)i/2, i%2) = std::atof (argv[i+9]); // loop parameters
    }

    BoxMatrix loopBox;
    AtomVector loopAtoms;
    Eigen::Matrix<float,3,3> loopMiller;
    loopMiller = Eigen::Matrix3f::Identity(); // 11-20 loop
    //loopMiller << 0.5, 0.0, -1.5,
    //  0.0, 1.0, 0.0,
    //  0.5, 0.0, 0.5;
    Eigen::Vector3f bLoop;
    bLoop << 1.0, 0.0, 0.0; // crystal coordinates
    printf ("... adding SIA loop with b=<x>\n");
    tmp = enki::create_sia_loop ( cell, Nx,
                                      loopMiller,
                                      Nloop,
                                      bLoop,
                                      loopAtoms, loopBox);
    printf ("+++ created SIA loop with %i self-interstitial atoms\n", tmp);
    printf ("    loop dimensions (box)\n");
    std::cout << loopBox << std::endl;
    atoms.insert(atoms.end(), loopAtoms.begin(), loopAtoms.end());

  }

  std::FILE * fp;
  fp = fopen ( "atoms.data", "w");
  enki::wrap_atoms_box(atoms, box);
  tmp = enki::write_lammps_data_file ( fp, atoms, box );
  printf ( "... done writing %i atoms to atoms.data\n", tmp );


  return 0;

}
