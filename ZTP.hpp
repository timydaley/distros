/****  
 
        Description: ...

  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Timothy Daley, Andrew D. Smith
 
  Authors: Philip J. Uren, Timothy Daley
 
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
  --------------------
 
  Known Bugs:    None
 
  Revision
  History:       None
         
  TODO:          None
 
****/

#ifndef ZTP_HPP
#define ZTP_HPP

#include <string>
#include <vector>
#include <exception>

/*****
 * @summary: ZTP exception class
 */
class ZTPError: public std::exception {
public:
  ZTPError (const char* msg){this->msg = msg;}
  virtual const char* what() const throw(){return this->msg;}
private:
  const char* msg;
};

/*****
 * @summary: ZTP class
 */
class ZeroTruncatedPoisson {
private:
  double lambda;
  static const double tolerance = 1e-20;
public:
  // constructors and destructors
  ZeroTruncatedPoisson(const double lambda);
  ZeroTruncatedPoisson();
       
  // mutators
  void estimateParams(const std::vector<double> &hist, bool debug = false);
       
  // inspectors
  double logLikelihood(const int obs) const;
  double pdf(const int y) const;
  double cdf(const int y) const;
  std::string toString() const;
  double get_lambda() const {return lambda;}
};

#endif
  
