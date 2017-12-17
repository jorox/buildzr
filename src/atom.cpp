#include "atom.h"
#include "Eigen/Dense"

Atom::Atom() {

  this->_id = -1;
  this->coords << 0., 0., 0.;
  this->ucc << 0, 0, 0;

}

Atom::Atom( int id,
            const Eigen::Vector3f& vscal,
            const Eigen::Vector3f& vreal){

  this->_id = id;
  this->ucc = vscal;
  this->coords = vreal;
}
