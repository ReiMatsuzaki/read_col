program main
  use Mod_UnitTest
  use Mod_BlockMat
  use Mod_BlockVec

  call run_test("SymVec", test_BlockVec)
  
  call run_test("BlockMat", test_BlockMat)
  call run_test("BlockMat_get", test_BlockMat_get)

  call run_test("block_matmul", test_block_matmul)
  
contains
  subroutine test_BlockVec()
    type(BlockVec) vec, vec2
    integer, parameter :: ifile = 13
    complex*16, allocatable :: vs(:)
    
    
    call BlockVec_new(vec, (/3, 2, 4/))
    call BlockVec_set_val(vec, 1, 1, (1.0d0, 0.0d0))
    call BlockVec_set_val(vec, 1, 2, (1.0d0, 1.0d0))
    call BlockVec_set_val(vec, 1, 3, (1.0d0, 2.0d0))

    ! call BlockVec_set_val(vec, 2, 1, (2.0d0, 0.0d0))
    ! call BlockVec_set_val(vec, 2, 2, (2.0d0, 1.0d0))
    call BlockVec_set_block(vec, 2, (/(2.0d0, 0.0d0), (2.0d0, 1.0d0)/))

    call BlockVec_set_val(vec, 3, 1, (3.0d0, 0.0d0))
    call BlockVec_set_val(vec, 3, 2, (3.0d0, 1.0d0))
    call BlockVec_set_val(vec, 3, 3, (3.0d0, 2.0d0))
    call BlockVec_set_val(vec, 3, 4, (3.0d0, 3.0d0))
    

    open(unit=ifile, file='test_symvec.dat')
    call BlockVec_write(vec, ifile)
    close(ifile)
    
    open(unit=ifile, file='test_symvec.dat')
    call BlockVec_new_read(vec2, ifile)
    close(ifile)

    call expect_eq("num_sym", 3, BlockVec_num_sym(vec2))
    call expect_eq("size(vec(1))", 3, BlockVec_block_size(vec2, 1))
    call expect_eq("vec(1)(2)", (1.0d0, 1.0d0), BlockVec_val(vec2, 1, 2))
    call Expect_Eq("vec2(2)(2)", (2.0d0, 1.0d0), BlockVec_val(vec2, 2, 2))
    call Expect_Eq("vec2(3)(3)", (3.0d0, 2.0d0), BlockVec_val(vec2, 3, 3))

    allocate(vs(BlockVec_block_size(vec2, 1)))
    call BlockVec_get_block(vec2, 1, vs)
    call expect_eq("vec(1)", (/(1.0d0,0.0d0), (1.0d0,1.0d0), (1.0d0,2.0d0)/), vs)
    deallocate(vs)

    allocate(vs(BlockVec_size(vec2)))
    call BlockVec_get_array(vec2, vs)
    call expect_eq("get_array", &
         (/(1.0d0, 0.0d0), &
         (1.0d0, 1.0d0), &
         (1.0d0, 2.0d0), &
         (2.0d0, 0.0d0), &
         (2.0d0, 1.0d0), &
         (3.0d0, 0.0d0), &
         (3.0d0, 1.0d0), &
         (3.0d0, 2.0d0), &
         (3.0d0, 3.0d0)/), vs)
    deallocate(vs)
    
  end subroutine test_BlockVec
  subroutine test_BlockMat()
    type(BlockMat) mat
    integer :: num_isym(4)    = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    integer num_i, num_j, num_ij(2)

    call BlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    call Expect_Eq("iblock_ijsym(1,1)", 1, mat % iblock_ijsym(1, 1))
    call Expect_Eq("iblock_ijsym(2,2)", 2, mat % iblock_ijsym(2, 2))
    call Expect_Eq("iblock_ijsym(1,2)", 5, mat % iblock_ijsym(1, 2))
    call Expect_Eq("offset(1)", 0, mat % offset_iblock(1))
    call Expect_Eq("idx_1", 1, BlockMat_index(mat, 1, 1, 1, 1))

    call BlockMat_set_val(mat, 1, 1, 1, 1, (1.0d0, 0.0d0))
    call BlockMat_set_val(mat, 1, 1, 1, 2, (1.0d0, 0.1d0))
    call BlockMat_set_val(mat, 1, 1, 2, 1, (1.0d0, 0.2d0))
    call BlockMat_set_val(mat, 1, 1, 2, 2, (1.0d0, 0.3d0))

    call BlockMat_set_val(mat, 2, 2, 1, 1, (2.0d0, 0.0d0))
    
    call BlockMat_set_val(mat, 3, 3, 1, 2, (3.0d0, 0.1d0))
    call BlockMat_set_val(mat, 3, 3, 2, 1, (3.0d0, 0.2d0))
    call BlockMat_set_val(mat, 3, 3, 2, 2, (3.0d0, 0.3d0))

    call BlockMat_set_val(mat, 1, 2, 1, 1, (5.0d0, 0.0d0))
    call BlockMat_set_val(mat, 1, 2, 1, 2, (5.0d0, 0.1d0))

    ! call BlockMat_show(mat)
    call Expect_Eq("1,1,1,1", (1.0d0, 0.0d0), BlockMat_val(mat, 1, 1, 1, 1))
    call Expect_Eq("1,1,1,2", (1.0d0, 0.1d0), BlockMat_val(mat, 1, 1, 1, 2))
    call Expect_Eq("1,2,1,2", (5.0d0, 0.1d0), BlockMat_val(mat, 1, 2, 1, 2))

    num_ij = BlockMat_block_size(mat, 1, 1)
    num_i = num_ij(1); num_j = num_ij(2)
    call expect_eq("num sym", 4, BlockMat_num_sym(mat))
    call expect_true("block size(1,1)", BlockMat_exist_block(mat, 1, 1))
    call Expect_Eq("block size(1,1)", 2, num_i)
    call Expect_Eq("block size(1,1)", 2, num_j)

    call expect_false("Block(1,3)", BlockMat_exist_block(mat, 1, 3))

    num_ij = BlockMat_block_size(mat, 1, 2)
    call expect_true("block size(1,2)", BlockMat_exist_block(mat, 1, 1))
    call Expect_Eq("block size(1,2)", 2, num_ij(1))
    call Expect_Eq("block size(1,2)", 1, num_ij(2))
    
    call BlockMat_delete(mat)
  end subroutine test_BlockMat
  subroutine test_BlockMat_get()
    type(BlockMat) mat
    integer :: num_isym(4) = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    complex*16 m11(2, 2), m12(2, 1)
    complex*16 a(10, 10), b(10, 10), x

    a = (0.0d0, 0.0d0)
    call BlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    x = (1.0d0,0.0d0); a(1,1) = x; call BlockMat_set_val(mat, 1, 1, 1, 1, x) 
    x = (1.0d0,0.1d0); a(1,2) = x; call BlockMat_set_val(mat, 1, 1, 1, 2, x) 
    x = (1.0d0,0.2d0); a(2,1) = x; call BlockMat_set_val(mat, 1, 1, 2, 1, x)
    x = (1.0d0,0.3d0); a(2,2) = x; call BlockMat_set_val(mat, 1, 1, 2, 2, x)
    
    x = (2.0d0,0.0d0); a(3,3) = x; call BlockMat_set_val(mat, 2, 2, 1, 1, x)
    
    x = (3.0d0,0.0d0); a(4,5) = x; call BlockMat_set_val(mat, 3, 3, 1, 2, x)
    x = (3.0d0,0.1d0); a(5,4) = x; call BlockMat_set_val(mat, 3, 3, 2, 1, x)
    x = (3.0d0,0.2d0); a(5,5) = x; call BlockMat_set_val(mat, 3, 3, 2, 2, x)

    x = (5.0d0,0.0d0); a(1,3) = x; call BlockMat_set_val(mat, 1, 2, 1, 1, x)
    x = (5.0d0,0.1d0); a(2,3) = x; call BlockMat_set_val(mat, 1, 2, 2, 1, x)

    call BlockMat_get_block(mat, 1, 1, m11)
    call Expect_Eq("m11_11", (1.0d0, 0.0d0), m11(1, 1))
    call Expect_Eq("m11_12", (1.0d0, 0.1d0), m11(1, 2))
    call Expect_Eq("m11_21", (1.0d0, 0.2d0), m11(2, 1))
    call Expect_Eq("m11_22", (1.0d0, 0.3d0), m11(2, 2))

    call BlockMat_get_block(mat, 1, 2, m12)
    call Expect_Eq("m12_11", (5.0d0, 0.0d0), m12(1, 1))
    call Expect_Eq("m12_21", (5.0d0, 0.1d0), m12(2, 1))
    
    call BlockMat_get_array(mat, b)
    call expect_eq("get_array", a, b)
    
  end subroutine test_BlockMat_get
  subroutine test_block_matmul()
    use Mod_BlockLinearAlgebra
    use Mod_Utils
    type(BlockVec) :: a, c
    type(BlockMat) :: B
    integer :: ist(3) = (/2, 1, 3/)
    complex*16 :: a1(2), c1(2)
    complex*16 :: a2(1), c2(1)
    complex*16 :: a3(3), c3(3)
    complex*16 :: b11(2, 2)
    complex*16 :: b22(1, 1)
    complex*16 :: b33(3, 3) 

    ! ---- compute vector a and matrix B ----
    call c_random(a1)
    call c_random(a2)
    call c_random(a3)
    call c_random(b11)
    call c_random(b22)
    call c_random(b33)    
    
    call BlockVec_new(a, ist)
    call BlockVec_set_block(a, 1, a1)
    call BlockVec_set_block(a, 2, a2)
    call BlockVec_set_block(a, 3, a3)
    
    call BlockMat_new(B, ist, (/1, 2, 3/), (/1, 2, 3/))
    call BlockMat_set_block(B, 1, 1, b11)
    call BlockMat_set_block(B, 2, 2, b22)
    call BlockMat_set_block(B, 3, 3, b33)

    ! ---- compute normal ----
    c = block_matmul(B, a)
    call BlockVec_get_block(c, 1, c1)
    call BlockVec_get_block(c, 2, c2)
    call BlockVec_get_block(c, 3, c3)

    call expect_eq("(Ba)(1)", c1, matmul(b11, a1))
    call expect_eq("(Ba)(2)", c2, matmul(b22, a2))
    call expect_eq("(Ba)(3)", c3, matmul(b33, a3))

    ! ---- compute transpose ----
    c = block_matmul(B, a, 't')
    call BlockVec_get_block(c, 1, c1)
    call BlockVec_get_block(c, 2, c2)
    call BlockVec_get_block(c, 3, c3)

    call expect_eq("(BT a)(1)", c1, matmul(transpose(b11), a1))
    call expect_eq("(BT a)(2)", c2, matmul(transpose(b22), a2))
    call expect_eq("(BT a)(3)", c3, matmul(transpose(b33), a3))
    
  end subroutine test_block_matmul
end program main
