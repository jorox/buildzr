#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"

struct UnitCell {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Eigen::Vector3f X, Y, Z;
  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> motif;
};

struct Atom {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  int _id;
  Eigen::Vector4i ucc; // unit-cell coordinates
  Eigen::Vector3f coords; // real coordinates

  Atom (int id): _id(id) {}
};

/**
   The function changes an atom index (int) into unit-cell corrdinates
 **/
void map_id_uc(int& idAtom, Eigen::Vector4i& uccAtom, Eigen::Matrix<int, 3,2>& boxUcLims, int nBasisAtoms){

  if (idAtom < 0){
    //printf("reverse mapping %i %i %i %i\n", uccAtom(0), uccAtom(1), uccAtom(2), uccAtom(3));
    //fflush(stdout);
    int Nx = boxUcLims(0,1) - boxUcLims(0,0);
    int Ny = boxUcLims(1,1) - boxUcLims(1,0);
    int Nz = boxUcLims(2,1) - boxUcLims(2,0);

    idAtom = uccAtom(0); //ib
    idAtom += (uccAtom(1) - boxUcLims(0,0)) * nBasisAtoms; // ix * Nb
    idAtom += (uccAtom(2) - boxUcLims(1,0)) * nBasisAtoms * Nx;
    idAtom += (uccAtom(3) - boxUcLims(2,0)) * nBasisAtoms * Nx * Ny;

  }
  else{
    int ib, ix, iy, iz; // indices in order of precedence
    int Nx, Ny, Nz; // number of unit cells
    int maxId; // max id possible

    Nx = boxUcLims(0,1) - boxUcLims(0,0); // number of unit cells aling X
    Ny = boxUcLims(1,1) - boxUcLims(1,0); // number of unit cells along Y
    Nz = boxUcLims(2,1) - boxUcLims(2,0); // number of unit cells aling Z

    maxId = Nx * Ny * Nz * nBasisAtoms; // max id possible

    ib = idAtom % nBasisAtoms; // fastest changing index
    ix = (int) idAtom / nBasisAtoms; // next fastest index
    ix = (ix % Nx) + boxUcLims(0,0);
    iy = (int) idAtom / (nBasisAtoms * Nx);
    iy = (iy % Ny) + boxUcLims(1,0);
    iz = (int) idAtom / (nBasisAtoms * Nx * Ny);
    iz = (iz % Nz) + boxUcLims(2,0);

    uccAtom << ib, ix, iy, iz;
  }
}

/**
   Use an atom's unit-cell coordinates to calculate its coordinates in real space
   The method accepts non-orthogonal unit-cells. However, it has not been tested
 **/
void map_uc_real(Atom& aatom, UnitCell& unit){
  aatom.coords << 0., 0., 0.;
  // basis coordinates
  aatom.coords += unit.motif[aatom.ucc(0)](0) * unit.X;
  aatom.coords += unit.motif[aatom.ucc(0)](1) * unit.Y;
  aatom.coords += unit.motif[aatom.ucc(0)](2) * unit.Z;
  // unit cell coordinates
  aatom.coords += aatom.ucc(1) * unit.X;
  aatom.coords += aatom.ucc(2) * unit.Y;
  aatom.coords += aatom.ucc(3) * unit.Z;
}



/**
   Write an atom list to a LAMMPS data file
   Only orthogonal box, 1 atom type, 1 header, and 1 orientation defined for the moment
 **/
void write_lammps_data_file( std::FILE* fstream, std::vector<Atom>& atomList, UnitCell& uc, Eigen::Matrix<float,3,3>& boxLims ){

  // count number of valid atoms
  int natoms = atomList.size();
  for (int i = 0; i < atomList.size(); ++i){
    if (atomList[i]._id == -1){ natoms--; }
  }

  // write the header
  fprintf ( fstream, "#ver .1, HCP-Zr [11-20], [0001], [-1100], metal units\n");
  fprintf ( fstream, "    %i atoms\n", natoms);
  fprintf ( fstream, "    1 atom types\n\n");
  fprintf ( fstream, "%1.6f %1.6f xlo xhi\n", boxLims(0,0), boxLims(0,1));
  fprintf ( fstream, "%1.6f %1.6f ylo yhi\n", boxLims(1,0), boxLims(1,1));
  fprintf ( fstream, "%1.6f %1.6f zlo zhi\n\n", boxLims(2,0), boxLims(2,1));
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

}

/**
   Create an edge dislocation with Burgers vector b = x, and n = z
   Atoms are not deleted to preserve fast id-ucc mapping.
   Instead their id is changed to -1 in which case they are not written out.
 **/
int create_edge_xz(std::vector<Atom>& atomList, UnitCell& uc, Eigen::Matrix<float, 3, 3>& boxLims, Eigen::Matrix<int, 3,2>& boxUcLims){

  //
  int Nx = boxUcLims(0,1) - boxUcLims(0,0);
  int Ny = boxUcLims(1,1) - boxUcLims(1,0);
  int Nz = boxUcLims(2,1) - boxUcLims(2,0);

  // find the mid-plane along the Z-direction
  int midz = boxUcLims(2,1) + boxUcLims(2,0);
  if (midz % 2 != 0){
    midz++;
  }
  midz = (int)midz / 2;

  // determine atoms to delete using reverse map_id_uc
  int nAtomsDelete = uc.motif.size() * Ny * (boxUcLims(2,1) - midz);
  int iadelete = -1; // index of atom to delete
  Eigen::Vector4i ucadelete; // unit-cell coordinates of atom to delete
  int ix = boxUcLims(0,1) - 1;

  // delete the atoms in at and above midz along the y-direction in the last unit-cell of X
  for (int iz = midz; iz < boxUcLims(2,1); ++iz){
    for (int iy = boxUcLims(1,0); iy < boxUcLims(1,1); ++iy){
      for (int ib = 0; ib < uc.motif.size(); ++ib){
        iadelete = -1;
        ucadelete << ib, ix, iy, iz;
        //printf ( "trying %i %i %i %i", ib, ix, iy, iz);
        map_id_uc( iadelete, ucadelete, boxUcLims, uc.motif.size());
        //printf ( " -- deleting %i %i %i %i\n", atomList[iadelete].ucc(0), atomList[iadelete].ucc(1), atomList[iadelete].ucc(2), atomList[iadelete].ucc(3) );
        atomList[iadelete]._id = -1;
      }
    }
  }

  Eigen::Vector4i ucatom;
  // strain the top and bottom parts
  for (int ia = 0; ia < atomList.size(); ++ia){
    if (atomList[ia].ucc(3) < midz) {
      atomsList[ia].coords(0)
    }




  return nAtomsDelete;

}


int main(int argc, char** argv){

  if (argc < 7){ //ensure that 6 number are given -- they will be cast as integers
    std::cout << "ERROR: must supply at least 6 integers to determine box limits" << std::endl;
    return 1;
  }

  /*
    lattice constants
   */
  double _ALAT = 3.23; // Angstrom
  double _CLAT = 5.165;

  /*
    Zirconium unit cell
    X = 1/3[11-20]
    Y = [0001]
    Z = [-1100]
    We have 4 basis atoms in this unit cell i.e. it is not the primitive; however, it is orthogonal
   */
  struct UnitCell zrstd;

  zrstd.X << _ALAT, 0.0, 0.0; // 1/3[11-20]
  zrstd.Y << 0.0, _CLAT, 0.0; // [0001]
  zrstd.Z << 0.0, 0.0, sqrt(3.0) * _ALAT; // [1-100]

  zrstd.motif.push_back(Eigen::Vector3f (0.0, 0.0, 0.0));   //A-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.0, 0.5, 1./3.)); //B-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.5, 0.0, 0.5));   //A-plane
  zrstd.motif.push_back(Eigen::Vector3f (0.5, 0.5, 0.5+1./3.)); //B-plane

  /*
    The box is defined by the repeat vectors Nx[:,i]
    The bottom-left and top-right corners are the multiple of the unit cell vectors with the number of unit cells
    along that direction. Note that in this way the bottom-left corner might NOT be the left-most point if the
    unit cell is not orthogonal for example.
   */
  Eigen::Matrix<int, 3, 2> Nx; // box dimensions in unit cell
  Eigen::Matrix<float, 3, 3> box_lim; // box limits i.e bottom-left and top-right corners
  int Natoms = 0;

  for (int i = 1; i<7; ++i){
    Nx((i - 1)/2, (i-1)%2) = atoi(argv[i]);
  }
  // calculate the blc and trc
  box_lim.block<3,1>(0,0) = Nx(0,0) * zrstd.X + Nx(1,0) * zrstd.Y + Nx(2,0) * zrstd.Z; //lower-left corner
  box_lim.block<3,1>(0,1) = Nx(0,1) * zrstd.X + Nx(1,1) * zrstd.Y + Nx(2,1) * zrstd.Z; //upper-right corner

  // calculate the number of atoms
  Natoms = (Nx(0,1) - Nx(0,0)) * (Nx(1,1) - Nx(1,0)) * (Nx(2,1) - Nx(2,0)) * zrstd.motif.size();

  // print information on the box
  printf ( "... box dimensions in lattice units: \n     X = [%i %i]\n     Y = [%i %i]\n     Z = [%i %i]\n",
           Nx(0,0), Nx(0,1), Nx(1,0), Nx(1,1), Nx(2,0), Nx(2,1));
  printf ( "    => %i atoms\n\n", Natoms);

  printf ( "... spacings along X, Y, Z: %1.5f, %1.5f, %1.5f\n", zrstd.X.norm(), zrstd.Y.norm(),zrstd.Z.norm());

  printf ( "... box dimensions in A: \n     X = %1.4f --> %1.4f\n     Y = %1.4f --> %1.4f\n     Z = %1.4f --> %1.4f\n\n",
           box_lim(0,0), box_lim(0,1), box_lim(1,0), box_lim(1,1), box_lim(2,0), box_lim(2,1));


  printf("... building crystal\n");
  std::vector<Atom> atoms;

  for (int i = 0; i < Natoms; ++i){
    atoms.push_back(Atom(i));
    map_id_uc ( i, atoms[i].ucc, Nx, zrstd.motif.size() ); // change atom unit-cell coordinates
    map_uc_real ( atoms[i], zrstd);
  }

  printf( "... done building crystal\n");
  if (argc > 7){
    int tmp = create_edge_xz ( atoms, zrstd, box_lim, Nx );
    printf ( "+++ created edge dislocation with b=x, n=z, %i atoms deleted \n", tmp);
  }


  std::FILE * fp;
  fp = fopen ( "zr.data", "w");
  write_lammps_data_file ( fp, atoms, zrstd, box_lim );
  printf ( "... done writing to zr.data\n" );







  /** //Test map_id_uc function
  printf ( "... testing <map_id_uc> function\n");
  int tNb = 3;
  Eigen::Matrix<int, 3, 2> testbox;
  Eigen::Vector4i tucAtom;

  testbox << 0, 3,
    -1, 1,
    0, 1;

  int tNtot = (testbox(0,1) - testbox(0,0)) * (testbox(1,1) - testbox(1,0)) * (testbox(2,1) - testbox(2,0)) * tNb;

  printf ("%i atoms:\n Nb = %i \n Nx = %i, %i\n Ny = %i, %i\n Nz= %i, %i", tNtot, tNb,
          testbox(0,0), testbox(0,1),
          testbox(1,0), testbox(1,1),
          testbox(2,0), testbox(2,1));

  for (int i = 0; i < tNtot; ++i){
    map_id_uc( i, tucAtom, testbox, tNb);
    printf( "\n%i --> %i %i %i %i", i, tucAtom(0), tucAtom(1), tucAtom(2), tucAtom(3));
  }
  std::cout << std::endl;
  **/

  return 0;

}
