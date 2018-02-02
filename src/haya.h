#ifndef HAYA_H
#define HAYA_H

#include "Eigen/Desne"
#include "unitcell.h"
#include <vector>
#include <string>
#include <fstream>

namespace haya{
  /*
    loads unit cell information from file
   */
  void load_unit_cell_from_file(std::FILE* fin,
                                UnitCell& cell);

  /*
    process a line and return an integer to determine what to do
  */
  int process_line (const std::string str,
                    const std::vector<std::string>& buildTypes,
                    std::vector<double>& params);

  /*
    split a string using delimeters and stores it in a vector
    WARNING: the contents of the vector will be destroyed
   */
  void split5 (const std::string& str,
               std::vector<std::string>& cont,
               const std::string& delims = " ");
}

#endif
