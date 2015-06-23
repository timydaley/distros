# FinitePoissonDirichletSampler:
#   a tool to sample species according to the two parameter
#   Poisson Dirichlet distribution
#   see Pitman & Yor, Annals of Prob, 1995 sec 9.1 or
#   Hansen & Pitman, Statistics & Probability Letters, 2000 eqs 13 & 14
#   
#    Copyright (C) 2015 Timothy Daley
# 
#    Authors: Timothy Daley
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.	If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from random import random
import argparse
import csv
import sys



def sampleFPDD(kappa, theta, n) : 
  counts = []
  counts.append(1)
  current_total_counts = 1
  while (current_total_counts < n) :
    # u is uniform rv in (0, 1)
    u = random()
    test_val = 0.0
    for i, count in enumerate(counts) :
      test_val += (count + kappa)/(current_total_counts + theta)
      if (u <  test_val) :
        counts[i] += 1
        break
    if (u > test_val) :
      counts.append(1)
    current_total_counts += 1
  return(counts)



def main():
  parser = argparse.ArgumentParser(description = 'Random species sample following the finite Poisson-Dirichlet process')
  parser.add_argument('-o', '--output_filename', dest = 'output_file',
                      type = argparse.FileType('w'), 
                      help = 'output file name')
  parser.add_argument('-k', '--kappa', dest = 'kappa',
                      default = 1.0, type = float,
                      help = 'kappa parameter, must be positive')
  parser.add_argument('-m', '--pop_size', dest = 'm',
                      default = 1000000, type = int,
                      help = 'size of population to be sampled from')
  parser.add_argument('-n', '--sample_size', dest = 'n',
                      default = 1000000, type = int,
                      help = 'number of individuals to sample')
  parser.add_argument('-s', '--seed', dest = 'seed',
                      default = None, 
                      help = 'seed for random number generator')
  parser.add_argument('-V', '--VERBOSE', action = "store_true",
                      dest = 'VERBOSE', default = False,
                      help = "run in verbose mode")
  args = parser.parse_args()

  theta = args.kappa*args.m

  if args.VERBOSE : 
    print("alpha = %s" % -args.kappa, file = sys.stderr)
    print("theta = %s" % theta, file = sys.stderr)
    print("n = %s" % args.n, file = sys.stderr)

  counts = sampleFPDD(args.kappa, theta, args.n)

  if args.VERBOSE :
    print("counts = ", file = sys.stderr)
    for count in counts : 
      print("%s" % count, file = sys.stderr)
    print("sum of counts = %s" % sum(counts), file = sys.stderr)

  assert(sum(counts) == args.n)
  for count in counts:
    args.output_file.write("%s\n" % count)
#end of main

if __name__ == '__main__':
  main()
    
