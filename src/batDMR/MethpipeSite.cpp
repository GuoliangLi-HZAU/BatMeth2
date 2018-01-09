/*
  Copyright (C) 2015 University of Southern California
  Authors: Andrew D. Smith

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with This program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "MethpipeSite.hpp"

#include <string>
#include <iostream>
#include <sstream>

//#include "smithlab_utils.hpp"

using std::string;

struct DNAmethException {
  DNAmethException(std::string m) : message(m) {}
  std::string what() const {return message;}
  std::string message;
};

std::istream &
operator>>(std::istream &in, MSite &s) {
  string line;
  if (!getline(in, line))
    return in;

  std::istringstream iss(line);
  if (!(iss >> s.chrom >> s.pos >> s.strand
        >> s.context >> s.meth >> s.n_reads))
    throw DNAmethException("bad methcounts file");
  return in;
}


string
MSite::tostring() const {
  std::ostringstream oss;
  oss << chrom << '\t'
      << pos << '\t'
      << strand << '\t'
      << context << '\t'
      << meth << '\t'
      << n_reads;
  return oss.str();
}
