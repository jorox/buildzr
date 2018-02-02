#include "haya.h"
#include "Eigen/Dense"
#include "unitcell.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

void split5(const std::string& str,
            std::vector<std::string>& cont,
            const std::string& delims = " ")
{
  cont.clear();
  boost::split(cont, str, boost::is_any_of(delims));
}

void haya::load_unit_cell_from_file(std::FILE* fin, UnitCell &cell)
{
  std::string line;
  Container words;
  int ic = -1;
  double sc = 1.0;

  while (std::getline(fin, line)){
    haya::split5(line, words);

    if (strcmp(words[0], "scale")==0) {sc = stof(words[1]); continue;}
    if (strcmp(words[0], "a1")==0){ic = 0;}
    if (strcmp(words[0], "a2")==0){ic = 1;}
    if (strcmp(words[0], "a3")==0){ic = 2;}
    if (strcmp(words[0], "basis")==0) {ic = 3;}

    if (ic == -1){
      std::cout << "ERROR: haya::load_unit_cell_from_file() undefined token: " << words[0] << std::endl;
    }
    if (ic > 0 and ic < 3){
        cell.ucv.block<3,1>(0,i) << std::stof(words[1]), std::stof(words[2]), std::stof(words[3]);
    }
    if (ic == 3){
      cell.basis.conservativeResize(Eigen::NoChange_t, cell.basis.cols()+1);
      cell.basis.block<3, 1>(0, cell.basis.cols()-1) << std::stof(words[1]),
        std::stof(words[2]), std::stof(words[3]);
    }
    ic = -1;
  }

  std::cout << "... loaded unit cell from file:" << std::endl;
  cell.print_me();
}

int haya::process_line(const std::string str,
                       const std::vector<std::string> &buildTypes,
                       std::vector<double> &params)
{

  if (str.empty()){ return -1;} // empty string
  std::vector<std::string> words;
  haya::split5(str, line);

  if ( std::strcmp(words[0], "build") != 0 ){ return -1;} //not an action line

  // look for the build-type
  for (int i = 0; i < buildTypes.size(); ++i){
    if (strcmp( words[1], buildTypes[i] ) == 0){
      params.clear();
      for (int j = 2; j < words.size(); ++j){
        params.push_back(stod(words[i]));
      }
      return i;
    }
  }

  std::cout << "WARNING: haya::process_line() unrecognized build command type " << words[1] << std::endl;
  return -1;
}
