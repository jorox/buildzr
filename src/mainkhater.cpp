/*!
  \file mainkhater.h
  \brief File containing the main function.

  The file creates a UnitCell, a 3x4 Matrix box, and a vector of atoms.
  The function uses haya namesapce for parsing and enki for creating crystals

  \var typedef Eigen::Matrix<double,3,4> BoxMatrix
  \brief A type definition for a 3D box.

  \var typedef std::vector<Atom> AtomVector
  \brief A type definition for a vector of atoms
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "atom.h"
#include "unitcell.h"
#include "crystal.h"
#include "enki.h"
#include "haya.h"

static constexpr double _PI_ = 3.1415926535897;
typedef Eigen::Matrix<double,3,4> BoxMatrix;
typedef std::vector<Atom> AtomVector;

/*!
 * \brief A function that joins two BoxMatrix along the Z-axis and concatenates two atom vectors
 * \param atomsM the first atoms vector
 * \param boxM the lower box (Mu)
 * \param atomsL the second atoms vector
 * \param boxL the upper box (Lambda)
 * \param atomsFull the joined atoms vector
 * \param boxFull the joined box
 * \todo
 *   Generalize the function to join along any direction
 *   Verify that the two boxes have the same origin
 *   Move the function to enki namespace
 */
void join (const AtomVector& atomsM,
             const BoxMatrix& boxM,
             const AtomVector& atomsL,
             const BoxMatrix& boxL,
             AtomVector& atomsFull,
             BoxMatrix& boxFull)
{
  // first join the atom lists
  atomsFull.reserve(atomsM.size()+atomsL.size());
  atomsFull.insert(atomsFull.end(), atomsM.begin(), atomsM.end());
  atomsFull.insert(atomsFull.end(), atomsL.begin(), atomsL.end());

  // join the boxes along the Z
  boxFull.block<3,2>(0,0) = 0.5 * (boxM.block<3,2>(0,0) + boxL.block<3,2>(0,0)); // The X and Y
  boxFull.block<3,1>(0,2) = boxM.block<3,1>(0,2) + boxL.block<3,1>(0,2); // the Z-axis
  boxFull.block<3,1>(0,3) = boxM.block<3,1>(0,3); // the origin = lower box
}
/*!
 * \brief Store 6-elements from a vector to a 3x2 matrix
 * This is used when reading the limits of a crystal from
 * \see haya::process_line()
 * \param vec is a vector that should contain 6 elements
 * \param mat the matrix to store the 6 values
 * \todo
 *   Make sure that the vector contains at least 6 elements
 */
void vec2mat (const std::vector<int> &vec, Eigen::Matrix<int, 3,2> &mat)
{
  mat << vec[0], vec[1], vec[2], vec[3], vec[4], vec[5];
}

int main(int argc, char** argv)
{

  if (argc < 2){ //ensure that 6 number are given -- they will be cast as integers
  std::cout << "ERROR: must supply an input file" << std::endl;
  return 1;
  }

  // program capabilities
  std::vector<std::string> xtalTypes;
  std::vector<int> params;
  xtalTypes.push_back("perfect");
  xtalTypes.push_back("edge");
  xtalTypes.push_back("screw");
  xtalTypes.push_back("sia");

  // main containers
  UnitCell cell0;
  AtomVector atoms;
  BoxMatrix box;

  std::ifstream ifs (argv[1], std::ifstream::in);

  // load unit cell from file
  haya::load_unit_cell_from_file(ifs, cell0);
  Eigen::Vector3d box_spc;
  cell0.size(box_spc);
  printf ( "... spacings (x,y,z): %1.5f, %1.5f, %1.5f\n",
           box_spc(0),
           box_spc(1),
           box_spc(2) );
  ifs.close();

  // edge and SIA boxes and atom containers
  Eigen::Matrix<int, 3, 2> Nx, Nloop; // box dimensions in unit cell
  int Natoms;
  int tmp;
  BoxMatrix boxMatMu, boxMatLambda, boxLoop;
  AtomVector atomsMu, atomsLambda, atomsLoop;
  Eigen::Matrix<double,3,3> hklLoop;
  Eigen::Vector3d bLoop;

  // open file to process commands
  ifs.open (argv[1], std::ifstream::in);
  std::string command;

  int ic;
  int nline = 0;
  while (std::getline(ifs, command)) {
    nline++;
    //std::cout << command << std::endl;
    ic = haya::process_line(command, xtalTypes, params);

    if (ic == 0){
        vec2mat( params, Nx);
        if (params.size() < 6){
          printf ("ERROR: in %s:%i need at least 6 integers to create perfect crystal\n",
                  argv[1], nline);
        }
        printf ( "... building %ix%ix%i perfect crystal\n",
                 Nx(0,1)-Nx(0,0), Nx(1,1)-Nx(1,0), Nx(2,1)-Nx(2,0));
        Natoms = enki::create_perfect(cell0, Nx, atoms, box);
        printf( "+++ crystal built successfully, %i atoms\n", Natoms);
        printf ( "    Box (Angstrom) [h1|h2|h2|origin]:\n");
        std::cout << box << std::endl;
        continue;
    }
    if (ic == 1){
      vec2mat( params, Nx);
      printf ("... building %ix%ix%i edge-dislocated crystal with b=y, n=z\n",
              Nx(0,1)-Nx(0,0), Nx(1,1)-Nx(1,0), Nx(2,1)-Nx(2,0));
      Natoms = enki::create_edge_xz (cell0, Nx, atomsMu, boxMatMu, atomsLambda, boxMatLambda );
      printf ( "+++ created edge dislocation with b=x, n=z, %i atoms \n", Natoms);
      printf ( "    Box-Mu (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsMu.size());
      std::cout << "    boxMat=\n" << boxMatMu << std::endl;
      printf ( "    Box-Lambda (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsLambda.size());
      std::cout << "    boxLambda=\n" << boxMatLambda << std::endl;
      join ( atomsMu, boxMatMu, atomsLambda, boxMatLambda, atoms, box );
      continue;
    }
    if (ic == 2){
      vec2mat( params, Nx);
      printf ("... building %ix%ix%i scrw-dislocated crystal with b=x, n=z\n",
              Nx(0,1)-Nx(0,0), Nx(1,1)-Nx(1,0), Nx(2,1)-Nx(2,0));
      Natoms = enki::create_screw_xz (cell0, Nx, atomsMu, boxMatMu, atomsLambda, boxMatLambda );
      printf ( "+++ created screw dislocation with b=x, n=z, %i atoms \n", Natoms);
      printf ( "    Box-Mu (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsMu.size());
      std::cout << "    boxMat=\n" << boxMatMu << std::endl;
      printf ( "    Box-Lambda (Angstrom) [h1|h2|h2|origin]: %i atoms\n", atomsLambda.size());
      std::cout << "    boxLambda=\n" << boxMatLambda << std::endl;
      join ( atomsMu, boxMatMu, atomsLambda, boxMatLambda, atoms, box );
      continue;
    }
    if (ic == 3){
      vec2mat(params, Nx);
      printf ("... adding SIA loop with b=<x>\n");
      haya::get_miller_sia( command, hklLoop);
      printf ("   Miller indices of the loop = \n");
      std::cout << hklLoop << std::endl;
      Eigen::Matrix<double,3,2> looptiles;
      looptiles.block<3,1>(0,0) = hklLoop.inverse() * Nx.block<3, 1>(0, 0).cast<double>();
      looptiles.block<3,1>(0,1) = looptiles.block<3,1>(0,0) + Nx.block<3,1>(0,1).cast<double>();
      std::cout << "*********DEBUG: looptiles = *************\n";
      std::cout << looptiles;
      std::cout << "\n*****************************************\n";
      Nloop = looptiles.cast<int>();
      bLoop << 1.0, 0.0, 0.0; // let burgers be equal to loop a1
      tmp = enki::create_sia_loop ( cell0,
                                    hklLoop,
                                    Nloop,
                                    bLoop,
                                    atomsLoop, boxLoop);
      printf ("+++ created SIA Miller loop with %i self-interstitial atoms\n", tmp);
      printf ("    loop dimensions (box)\n");
      std::cout << boxLoop << std::endl;;
      atoms.insert(atoms.end(), atomsLoop.begin(), atomsLoop.end());
      continue;
    }

  }

  std::FILE * fp;
  fp = fopen ( "atoms.data", "w");
  enki::wrap_atoms_box(atoms, box);
  tmp = enki::write_lammps_data_file ( fp, atoms, box );
  printf ( "... done writing %i atoms to atoms.data\n", tmp );

  return 0;

}
