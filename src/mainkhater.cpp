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
    printf ( "+++  %i atoms\n", Natoms);

    std::FILE * fp;
    fp = fopen ( "zr.data", "w");
    int tmp = write_lammps_data_file ( fp, atoms, boxMat );
    printf ( "... done writing %i atoms to zr.data\n", tmp );
  }


  if (argc > 7){
    Eigen::Matrix<float,3,4> boxMat1;
    std::vector<Atom> atoms1;
    int tmp;
    int dcase = -1;
    if (strcmp(argv[7],"e") == 0){
      tmp = create_edge_xz (zrstd_uc, Nx, atoms, boxMat, atoms1, boxMat1 );
      printf ( "+++ created edge dislocation with b=x, n=z, %i atoms \n", tmp);
    dcase = 1;
    }
    else{
      if ( strcmp(argv[7], "s") == 0 ){
        tmp = create_screw_xz(zrstd_uc, Nx, atoms, boxMat, atoms1, boxMat1);
        std::cout << "ddd boxMat=\n" << boxMat << std::endl;
        std::cout << "ddd boxMat1=\n" << boxMat1 << std::endl;
        printf ( "+++ created screw dislocation with b=x, n=z,\n");
        printf ( "    %i atoms:\n", tmp);
        dcase = 2;
      }

      else{
        if (strcmp(argv[7], "l") == 0 ) {
          printf ( "... building perfect crystal\n");
          int Natoms = create_perfect(zrstd_uc, Nx, atoms, boxMat);
          printf( "    crystal built successfully\n");
          printf ( "+++ Box (Angstrom) [h1|h2|h2|origin]:\n");
          std::cout << boxMat << std::endl;
          printf ( "+++  %i atoms\n", Natoms);


          Eigen::Matrix<float,3,2> loopTiles;
          for (int i = 0; i < 6; ++i){ loopTiles(i/2, i%2) = atof(argv[i+8]); }
          Eigen::Matrix<float,3,3> loopMiller;
          //loopMiller = Eigen::Matrix3f::Identity(); // 11-20 loop
          loopMiller << 0.5, 0.0, -1.5,
            0.0, 1.0, 0.0,
            0.5, 0.0, 0.5;
          Eigen::Vector3f bLoop;
          bLoop << 1.0, 0.0, 0.0; // crystal coordinates
          tmp = create_sia_loop ( zrstd_uc, Nx,
                                  loopMiller,
                                  loopTiles,
                                  bLoop,
                                  atoms1, boxMat1);
          printf ("... created SIA loop with %i atoms\n", tmp);
          dcase = 3;

        }
      }
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
    if (dcase < 3){
      boxFull.block<3,2>(0,0) = 0.5 * (boxMat.block<3,2>(0,0) + boxMat1.block<3,2>(0,0));
      boxFull.block<3,1>(0,2) = boxMat.block<3,1>(0,2) + boxMat1.block<3,1>(0,2);
      boxFull.block<3,1>(0,3) = boxMat.block<3,1>(0,3);
    }
    else{
      boxFull = boxMat;
    }
    enki::wrap_atoms_box(atomsFull, boxFull);
    tmp = write_lammps_data_file (fpp, atomsFull, boxFull );
    printf ( "    done writing %i atoms to zr.data\n", tmp);
  }

  return 0;

}
