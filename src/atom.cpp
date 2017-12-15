#include "atom.h"
#include "Eigen/Dense"

Atom::Atom() {

  this->_id = -1;
  this->coords << 0., 0., 0.;
  this->ucc << 0, 0, 0, 0;

}

Atom::Atom( int id, Eigen::Vector4i& vfrac, Eigen::Vector3f& vreal){
  this->_id = id;
  this->ucc = vfrac;
  this->coords = vreal;
}
