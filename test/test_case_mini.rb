
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

class TC_small_Read_ccolumbus < Test::Unit::TestCase
  
  def setup

    # obj build from small matrix
    @obj_mini = Read_ccolumbus.new ["S", "P", "D"], [["P", "S"], ["D", "D"]]
    s_basis = [{zeta: Complex(1.0, -2.0), linear_comb: {"XX" => 1.0}},
               {zeta: Complex(2.0, -3.0), linear_comb: {"XX" => 1.0}}]
    p_basis = [{zeta: Complex(1.0, -2.0), linear_comb: {"XX" => 1.0}},
               {zeta: Complex(2.0, -3.0), linear_comb: {"XX" => 1.0}},
               {zeta: Complex(2.0, -1.0), linear_comb: {"XX" => 1.0}}]
    d_basis = [{zeta: Complex(1.0, -2.0), linear_comb: {"XX" => 1.0}},
               {zeta: Complex(2.0, -3.0), linear_comb: {"XX" => 1.0}}]
    @obj_mini.symmetry_adapeted_basis_sym_basis = 
      {"S"=>s_basis, "P"=>p_basis, "D"=>d_basis}
    @obj_mini.configuration_I = [
      {:mo_symmetry=>"P", :index=>1}, {:mo_symmetry=>"S", :index=>1},
      {:mo_symmetry=>"P", :index=>3}, {:mo_symmetry=>"S", :index=>1},
      {:mo_symmetry=>"D", :index=>1}, {:mo_symmetry=>"P", :index=>2}]
    @obj_mini.mo_sym_ni = {
      "S" => NArray.to_na([[0.1, 0.2], [0.2, 0.5]]),
      "P" => NArray.to_na([[0.1, 0.2, 0.3], [0.2, -0.2, 0.5], [0.3, -0.5, 1.2]]),
      "D" => NArray.to_na([[0.1, 0.2], [0.2, 0.5]])}
    @obj_mini.ci_symsym_Iij = {
      ["P", "S"] => NArray.to_na(
          [[[1.2, 1.3, 2.4], [1.3, 0.4, 0.5], [1.31, 0.43, 0.52],],
           [[1.2, 1.3, -2.4],[1.1, 0.2,-0.5], [-1.3, -0.4, -0.15],]]),
      ["D", "P"] => NArray.to_na([
           [[1.2, 1.3, 1.3] ,[0.4, 2.3, 0.1]],
           [[1.2, 1.3, 1.3], [0.4, 2.3, 0.1]],
           [[3.2,-2.3, 6.6], [-3.4, 0.3, -2.1]]])
    }
    
    @obj_mini.delta_ene_I = NArray.to_na([0.5, 1.65, 2.01])
    @obj_mini.l_tdm_I = NArray.to_na([-2.2, 1.65, 0.01])
    @obj_mini.v_tdm_I = NArray.to_na([-1.2, 1.25, 0.11])    
  end

  def test_first_order
    fo_I = @obj_mini.compute_first_order_coef_length(1.5)
    assert_complex_in_epsilon(fo_I[1], Complex(11.0, 0.0))        
  end
  def test_c_Iij_fo_I_ij
    c_Iij_fo_I_ij = @obj_mini.compute_c_Iij_fo_I_ij(1.5, "P")
    assert_complex_in_epsilon(c_Iij_fo_I_ij[0,1], Complex(16.8929, 0.0), 0.001)            
  end
  def test_no_matrix_based_on_mo 
    no_ij = @obj_mini.compute_no_matrix_based_mo(1.5, "P")
    assert_complex_in_epsilon(no_ij[0,1], Complex(201.372, 0.0), 0.001)        
  end
  def test_no_orbital
    @obj_mini.test_show_mini()
    n_i, uinv_ij, coef_ni = @obj_mini.compute_no_orbital(1.5, "P")
    assert_complex_in_epsilon(n_i[2, 2], Complex(645.317, 0.0), 0.01)
    assert_complex_in_epsilon(uinv_ij[1,2], Complex(-0.333245, 0.0), 0.01)
    assert_complex_in_epsilon(coef_ni[0,1], Complex(-0.317271, 0.0), 0.01)
#    assert_complex_in_epsilon(coef_ni[1,0], Complex(-0.317271, 0.0), 0.01)        
  end

end

