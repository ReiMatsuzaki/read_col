#!/usr/bin/ruby

ccol = File.expand_path "~/ccolumbus/bin/"
work = "/scr/matsuzak/he_mini/"
out  = Dir::pwd + "/out/"

Dir::mkdir out  if not File.exist? out
Dir::mkdir work if not File.exist? work

system "cp *.in #{work}"
Dir::chdir work

puts "integrals..."
system "#{ccol}/molint < int.in > #{out}/int.out"
puts "dipoles..."
system "#{ccol}/cdip < int.in > cdip.out"
puts "pkfl..."
system "#{ccol}/pkfl < dum.in > Pkfl.out"
puts "scf..."
system "#{ccol}/cscf < scf.in > #{out}/scf.out"
puts "ci 1..."
system "#{ccol}/molci < ci1.in > #{out}/ci1.out"
system "mv CSF  CSFG"
system "mv CIVEC CIVECG"
puts "ci 2..."
system "#{ccol}/molci < ci2.in > #{out}/ci2.out"

system "bin_civec #{work}/CIVEC    > #{out}/civec.out"
system "bin_civec #{work}/CIVECG 1 > #{out}/civecG.out"
system "cp        #{work}/phoxsec #{out}/"
system "cp        #{work}/MOCOEF #{out}/"

