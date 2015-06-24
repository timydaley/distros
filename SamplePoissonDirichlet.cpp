/*    SamplePoissonDirichlet
 *
 *    Copyright (C) 2015 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <queue>
#include <string>
#include <sys/types.h>
#include <fstream>
#include <iostream>
#include <sstream>



#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <MappedRead.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>

#include "PoissonDirichlet.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;

    bool VERBOSE = false;
    double kappa = 0.5;
    size_t pop_size = 100000;
    size_t sample_size = 100000;



    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("SamplePoissonDirichlet", 
			   "sample species counts that follow a finite Poisson Dirichlet process");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
    		      false , outfile);
    opt_parse.add_opt("kappa",'k', "kappa, default = 0.5", false, kappa);
    opt_parse.add_opt("pop_size",'m',"population size to sample from, default = 100,000", 
		      false, pop_size);
    opt_parse.add_opt("sample_size", 'n', "number of individuals to capture, default = 100,000", 
		      false, sample_size);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    /**********************************************************************/
    
    const double sigma = -kappa;
    const double theta = kappa*pop_size;
    PDD real_distro(theta, sigma);
    vector<double> sample_counts;
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); 
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    real_distro.sample_PoissDir_counts(rng, sample_size, sample_counts);
    cerr << "number of counts = " << sample_counts.size() << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    for(size_t i = 0; i < sample_counts.size(); i++)
      out << sample_counts[i] << endl;

  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
