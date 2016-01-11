#!/usr/bin/ruby

require '/home/matsuzak/ccolumbus/script/make_input.rb'
require "complex"

InputFiles do
  int_input "int.in"
  scf_input "scf.in"
  dum_input "dum.in"
  ci_input  "ci1.in", "ci2.in"

  atom "He", 2.0, [[0.0, 0.0, 0.0]]
  
  symmetry  "atom"
  
  basis_set "He", "He/mini", {1=>1,2=>2,3=>3,4=>4}
  
  scf      "scfsd"
  w_range  "59.8000", "61.9000"
  ci       "cisd1", "cisd2"
end



