#include "haya.h"
#include "Eigen/Dense"
#include "unitcell.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

void haya::split5(const std::string& str,
            std::vector<std::string>& cont,
            const std::string& delims)
{
  cont.clear();
  boost::split(cont, str, boost::is_any_of(delims));
}

void haya::load_unit_cell_from_file(std::ifstream &fin, UnitCell &cell)
{
  std::string line;
  std::vector<std::string> words;
  int ic = -1;
  double sc = 1.0;

  while (std::getline(fin, line)){
    haya::split5(line, words);

    if ( words[0] == "#" ) { continue;}
    if ( words[0] == "scale" ) {sc = stod(words[1]); continue;}
    if ( words[0] == "name" ) {cell.name = words[1]; continue;}
    if ( words[0] == "a1" ){ic = 0;}
    if ( words[0] == "a2" ){ic = 1;}
    if ( words[0] == "a3" ){ic = 2;}
    if ( words[0] == "basis" ) {ic = 3;}

    if (ic == -1){
      std::cout << "ERROR: haya::load_unit_cell_from_file() undefined token: " << words[0] << std::endl;
    }
    if (ic > -1 and ic < 3){
        cell.ucv.block<3,1>(0,ic) << sc * std::stod(words[1]),
          sc * std::stod(words[2]), sc * std::stod(words[3]);
    }
    if (ic == 3){
      cell.basis.conservativeResize(3, cell.basis.cols()+1);
      cell.basis.block<3, 1>(0, cell.basis.cols()-1) << std::stod(words[1]),
        std::stod(words[2]), std::stod(words[3]); // the fractional coords
      cell.types.push_back(std::stoi(words[4])); // the type
    }
    ic = -1;
  }

  std::cout << "... loaded unit cell from file:" << std::endl;
  cell.print_me();
}

int haya::process_line(const std::string str,
                       const std::vector<std::string> &buildTypes,
                       std::vector<int> &params)
{
  if (str.empty()){ return -1;} // empty string
  std::vector<std::string> words;
  haya::split5(str, words);

  if ( words[0] != "build" and words[0] != "transform"){ return -1;} //not an action line
  if ( words[0] == "build" and words.size() < 8 ) {
    std::cout << "ERROR: Need at least 6 integers for build command\n";
    std::cout << str;
  }
  if ( words[0] == "transform" and words.size() < 4) {
    printf ("ERROR: Need at least 3 miller indices for a transform command\n");
    printf ( str.c_str() );
  }
  else{
    if (words[0] == "transform") {return 99;}
  }
  // look for the build-type
  params.clear();
  for (int i = 0; i < buildTypes.size(); ++i){
    if (  words[1] == buildTypes[i]){
      for (int j = 2; j < 8; ++j){params.push_back(std::stoi(words[j]));} //this wont work for a transform command
      return i;
    }
  }

  std::cout << "WARNING: haya::process_line() unrecognized build command type ";
  std::cout << words[1] << std::endl;
  return -1;
}

void haya::get_miller_sia( const std::string command,
                           Eigen::Matrix3d& hkls)
{
  std::vector<std::string> words;
  haya::split5(command, words);
  if (words.size() < 11){
    hkls << Eigen::Matrix3d::Identity();
  }
  else{
    for (int i=8; i < 11; ++i){
      hkls.block<3,1>(0,i-8) << str2miller ( words[i] );
    }
  }
}

void haya::get_miller_transform( const std::string command,
                                 Eigen::Matrix3d& hkls){
  std::vector<std::string> words;
  haya::split5(command, words);
  if (words.size() < 4){
    hkls << Eigen::Matrix3d::Identity();
  }
  else{
    for (int i=1; i < 4; ++i){
      hkls.block<3,1>(0,i-1) << str2miller ( words[i] );
    }
  }
}

Eigen::Vector3d haya::str2miller(const std::string str){
  std::string millerStr (str);
  boost::trim (millerStr);

  // find the braces
  std::size_t openBrace = millerStr.find("[");
  std::size_t closeBrace = millerStr.find( "]" );
  //std::cout << "position of '[' = " << openBrace << std::endl;
  //std::cout << "position of ']' = " << closeBrace << std::endl;

  // get the scale
  double scale = 1.0;
  if (openBrace != std::string::npos and openBrace >0){
    scale = std::stod ( millerStr.substr(0, openBrace) );
  }
  //std::cout << "scale = " << scale << std::endl;

  // get the indices
  std::string indx (millerStr.substr(openBrace!=std::string::npos? openBrace+1 : 0,
                                     closeBrace!=std::string::npos? closeBrace-openBrace-1:millerStr.length()));
  //std::cout << "indices = " << indx << std::endl;

  int n = indx.length();
  //if (n < 2 or n > 6){ }

  char letters[n+1];
  std::strcpy(letters, indx.c_str());

  // change to a vector
  Eigen::Vector3d hkl(1.,1.,1.);
  int imiller = 0;

  for (int i=0; i<n; ++i){
    if (isdigit(letters[i])){
      hkl(imiller) *= double( (int)letters[i] - 48 );
      imiller++;
    }
    if (letters[i] == 'm'){
      hkl(imiller) = -1.0;
    }
  }
  //std::cout << "hkl = ";
  hkl *= scale;

  return hkl;
}
