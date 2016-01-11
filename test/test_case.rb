
require "~/src/prog/cCOLUMBUS/script/read_ccol2/read_ccol2.rb"
require 'test/unit'
require "complex"
require "narray"

# Unit test script for read_ccol2.rb.
# Read the output of computation result of Helium using mini basis set.

class Test::Unit::TestCase
  def assert_complex_in_epsilon(actual, expected, epsilon=0.001, message="")
    assert_msg = message + "\n" + 
      "actual: #{actual}\n" +
      "expected: #{expected}\n" +
      "Error: #{(expected-actual).abs}\n" + 
      "Eps: #{epsilon}"
    assert((actual - expected).abs < epsilon, assert_msg)
  end
end

class TC_Read_ccolumbus < Test::Unit::TestCase
  
  def setup

    # obj build from reading ccolumbus file.
    sym_array =  ["S", "Pz", "Px", "Py", "Dxy", "Dxz", "Dyz", "xyz"]
    symsym_array = [["Pz", "S"], ["Dxz", "Px"], ["Dyz", "Py"],["xyz", "Dxy"]]
    @obj = Read_ccolumbus.new sym_array, symsym_array
    @obj.read_int("out/int.out")
    @obj.read_mo("out/MOCOEF")
    @obj.read_configurations("out/ci2.out")
    @obj.read_ci "out/civec.out"
    @obj.read_photo_excitation("out/phoxsec")

  end

  # test for reader
  def test_orbital_symmetry_array
    assert(@obj.orbital_symmetry_array != nil, "nil?")    
    assert(@obj.orbital_symmetry_array[1] == "Pz", "Second element is not Pz")
  end
  def test_symmetry_adapeted_basis_sym_basis
    assert(@obj.symmetry_adapeted_basis_sym_basis != nil, "nil?")
    assert(@obj.symmetry_adapeted_basis_sym_basis["Pz"]  != nil, "Pz")

    calc = @obj.symmetry_adapeted_basis_sym_basis["Pz"][1][:zeta]
    exact=Complex(1.5, -0.5)
    assert_complex_in_epsilon(calc, exact)
  end
  def test_mo_sym_ni
    assert_complex_in_epsilon(@obj.mo_sym_ni["Pz"][1,3],
                              Complex(3.439096, 1.797644))

  end
  def test_configuration_I
    assert_equal(@obj.configuration_I[11],
                 [{:mo_symmetry=>"Pz", :index=>1}, {:mo_symmetry=>"S", :index=>7}])
  end
  def test_ci_symsym_Iij
    actual = @obj.ci_symsym_Iij[["Pz", "S"]][3, 2, 0]
    expected  = Complex(1.473320881503251 * 0.01, 3.955113788243111 * 0.001)
    
    assert_complex_in_epsilon(actual, expected, 0.0001, "4 th excited state. 1S, 3Pz")
  end
  def test_excitation
    assert_complex_in_epsilon(@obj.delta_ene_I[3], Complex(2.22546, -0.01092))
    assert_complex_in_epsilon(@obj.l_tdm_I[3], Complex(0.00973592, 0.00193766))
    assert_complex_in_epsilon(@obj.v_tdm_I[3], Complex(0.05243865, 0.00967204))    
  end

end

