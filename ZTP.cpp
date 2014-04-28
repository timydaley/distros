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

#include "ZTP.hpp"
#include <vector>
#include <numeric>
#include <sstream>
#include <math.h>
#include <iostream>
#include <limits>
#include <gsl/gsl_sf.h>

using std::stringstream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

/*****
 * @summary: ZTP constructors
 */
ZeroTruncatedPoisson::ZeroTruncatedPoisson(const double lambda) : lambda(lambda) {}
ZeroTruncatedPoisson::ZeroTruncatedPoisson() : lambda(-1) {}
       
/*****
 * @summary: learn ZTP parameters
 */
static inline double
movement(const double a, const double b) {
  return fabs((a - b)/std::max(a, b)); //delta
}

static inline double
lambda_score_funct(const double mean, const double lambda){
  return(lambda + mean*exp(-lambda) - mean);
}

static double
trunc_estim_param_bisec(const double mean, const double tol){

  double lambda_low = mean-1;
  double lambda_high = mean;
  double lambda_mid = mean - 0.5;
  
  double diff = std::numeric_limits<double>::max();
  double prev_val = std::numeric_limits<double>::max();
  
  while(movement(lambda_high, lambda_low) > tol &&
         diff > tol){
    
    lambda_mid = (lambda_low + lambda_high)/2;
    
    const double mid_val = lambda_score_funct(mean, lambda_mid);
    
    if (mid_val < 0) lambda_low = lambda_mid;
    else lambda_high = lambda_mid;
    
    diff = fabs((prev_val - mid_val)/std::max(mid_val, prev_val));
    
    prev_val = mid_val;
  }
  return lambda_mid;
}

void ZeroTruncatedPoisson::estimateParams(const vector<double> &hist, bool debug) {
  double mu = 0.0;
  for(size_t i = 0; i < hist.size(); i++)
    mu += i*hist[i];
  mu = mu/accumulate(hist.begin(), hist.end(), 0.0);
  this->lambda = 
  trunc_estim_param_bisec(mu, this->tolerance);
}
       
/*****
 * @summary: ZTP prob. mass function
 * @note: probably has overflow issues with large y due to gsl_sf_fact(y)
 */
double ZeroTruncatedPoisson::pdf(const int y) const {
  if (y==0) return 0;
  return (exp(-this->lambda) * pow(this->lambda, y) /
	  gsl_sf_fact(y)) / (1 - exp(-this->lambda));
}

/*****
 * @summary: ZTP prob. mass function in log space
 */
double ZeroTruncatedPoisson::logLikelihood(const int y) const {
  if (y==0) return -INFINITY;
  return -lambda + (y*log(lambda)) 
	  - gsl_sf_lnfact(y) - log(1-exp(-lambda));
}

/*****
 * @summary: ZTP cummulative mass function
 */
double ZeroTruncatedPoisson::cdf(const int y) const {
  double total = 0;
  for (int i=0; i<=y; i++) {
    total += exp(this->logLikelihood(i));
  }
  return total;
}

/*****
 * @summary: get string representation of this ZTP object
 */
string ZeroTruncatedPoisson::toString() const {
  stringstream ss;
  ss << "lambda: " << this->lambda;
  return ss.str();
}
  
