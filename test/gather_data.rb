
require "~/src/prog/cCOLUMBUS/script/read_ccol2/read_ccol2.rb"

obj = Read_ccolumbus.new ["S", "Px", "Py", "Pz", "Dxy", "Dxz", "Dyz", "xyz"],  [["Pz", "S"], ["Dxz", "Px"], ["Dyz", "Py"],["xyz", "Dxy"]]

puts "reading files"
puts ""

obj.read_int("out/int.out")
obj.read_mo("out/MOCOEF")
obj.read_configurations("out/ci2.out")
#obj.configuration_array.each { |conf| p conf }
obj.read_ci "out/civec.out"
obj.read_photo_excitation("out/phoxsec")



puts ""
puts "====== test_show ========="
puts ""
obj.test_show

p obj.orbital_symmetry_array

omat = obj.compute_natural_orbital_matrix_based_mo(2.5, ["Pz", "S"])
p omat.shape
p omat[0,0]
