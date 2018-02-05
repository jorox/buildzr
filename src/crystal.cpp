#include "crystal.h"
#include "unitcell.h"
#include <vector>
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <iostream>
#include <stdio.h>

Crystal::Crystal(){
  this->_id = -1;
  this->uc = 0;
  this->tiles << 0,0,0,0,0,0;
}

Crystal::Crystal(int cid,
                 const Eigen::Matrix<int,3,2>& ctiles,
                 UnitCell* cuc){
  this->_id = cid;
  this->uc = cuc;
  this->tiles = ctiles;
}

int Crystal::build(std::vector<Atom >& atomList) const{
  int natoms = uc->generate(this->tiles, atomList);
  return natoms;
}

int Crystal::number_of_atoms() const{
  return ( ( this->tiles(0,1) - this->tiles(0,0) ) *
           ( this->tiles(1,1) - this->tiles(1,0) ) *
           ( this->tiles(2,1) - this->tiles(2,0) ) *
           this->uc->number_of_atoms()
           );
}

void Crystal::print_me() const{

  std::cout << "Crystal: %i" << std::endl;
  std::cout << "motif (cols) in A =\n" << this->uc->ucv << std::endl;
  std::cout << "basis (cols) fractional = \n" << this->uc->basis << std::endl;
}
