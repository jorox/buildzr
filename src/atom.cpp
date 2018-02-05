#include "atom.h"
#include "Eigen/Dense"

Atom::Atom() {

  this->_id = -1;
  this->coords << 0., 0., 0.;
  this->ucc << 0, 0, 0;

}

Atom::Atom( int id,
            int typ,
            const Eigen::Vector3d& vscal,
            const Eigen::Vector3d& vreal){

  this->_type = typ;
  this->_id = id;
  this->ucc = vscal;
  this->coords = vreal;
}
