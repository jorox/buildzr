#include "eigen/dense"

class UnitCell{

 public:
  /**
     Construct a crystalline system using 3 vectors and a motif
     The 3 lattice vectors should be specified in the standard e1e2e3 space. Units in Angstrom. They can be orthogonal or
     non-orthogonal. Care should be taken in case that the system is not orthonormal
   **/
  UnitCell (Eigen::Vector3f& v1, Eigen::Vector3f& v2, Eigen::Vector3f& v3, Eigen::Matrix<float, 3, d> motif);
  /**
     Return the volume of the unit cell
   **/
  double volume();
  /**
     Return the length of the three unit cell vectors
   **/
  void size(Eigen::Vector3f& );


 private:

  Eigen::Vector3f _a1;
  Eigen::Vector3f _a2;
  Eigen::Vector3f _a3;
  Eigen::Matrix<float, 3, d> _motif


}
