program main
  use Mod_AoInts
  use Mod_UnitTest  
  implicit none

  call run_test("BlockMat", test_BlockMat)
  call run_test("get_as_array", test_BlockMat_get_as_array)
  ! call run_test("AOINTS", test_AoInts)
  call run_test("DIPINTS", test_DIPINTS)
  
contains
  subroutine test_BlockMat()
    type(BlockMat) mat
    integer :: num_isym(4)    = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    integer num_i, num_j, num_ij(2)
    logical have_val

    call BlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    call Expect_Eq("iblock_ijsym(1,1)", 1, mat % iblock_ijsym(1, 1))
    call Expect_Eq("iblock_ijsym(2,2)", 2, mat % iblock_ijsym(2, 2))
    call Expect_Eq("iblock_ijsym(1,2)", 5, mat % iblock_ijsym(1, 2))
    call Expect_Eq("offset(1)", 0, mat % offset_iblock(1))
    call Expect_Eq("idx_1", 1, BlockMat_index(mat, 1, 1, 1, 1))

    call BlockMat_set_element(mat, 1, 1, 1, 1, (1.0d0, 0.0d0))
    call BlockMat_set_element(mat, 1, 1, 1, 2, (1.0d0, 0.1d0))
    call BlockMat_set_element(mat, 1, 1, 2, 1, (1.0d0, 0.2d0))
    call BlockMat_set_element(mat, 1, 1, 2, 2, (1.0d0, 0.3d0))

    call BlockMat_set_element(mat, 2, 2, 1, 1, (2.0d0, 0.0d0))
    
    call BlockMat_set_element(mat, 3, 3, 1, 2, (3.0d0, 0.1d0))
    call BlockMat_set_element(mat, 3, 3, 2, 1, (3.0d0, 0.2d0))
    call BlockMat_set_element(mat, 3, 3, 2, 2, (3.0d0, 0.3d0))

    call BlockMat_set_element(mat, 1, 2, 1, 1, (5.0d0, 0.0d0))
    call BlockMat_set_element(mat, 1, 2, 1, 2, (5.0d0, 0.1d0))

    ! call BlockMat_show(mat)
    call Expect_Eq("1,1,1,1", (1.0d0, 0.0d0), BlockMat_element(mat, 1, 1, 1, 1))
    call Expect_Eq("1,1,1,2", (1.0d0, 0.1d0), BlockMat_element(mat, 1, 1, 1, 2))
    call Expect_Eq("1,2,1,2", (5.0d0, 0.1d0), BlockMat_element(mat, 1, 2, 1, 2))

    num_ij = BlockMat_block_size(mat, 1, 1)
    num_i = num_ij(1); num_j = num_ij(2)
    call expect_true("block size(1,1)", have_val)
    call Expect_Eq("block size(1,1)", 2, num_i)
    call Expect_Eq("block size(1,1)", 2, num_j)

    call expect_false("Block(1,3)", BlockMat_exist_block(mat, 1, 3))

    num_ij = BlockMat_block_size(mat, 1, 2)
    call expect_true("block size(1,2)", have_val)
    call Expect_Eq("block size(1,2)", 2, num_ij(1))
    call Expect_Eq("block size(1,2)", 1, num_ij(2))
    
    call BlockMat_delete(mat)
  end subroutine test_BlockMat
  subroutine test_BlockMat_get_as_array()
    type(BlockMat) mat
    integer :: num_isym(4) = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    complex*16 m11(2, 2), m12(2, 1)

    call BlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    call BlockMat_set_element(mat, 1, 1, 1, 1, (1.0d0, 0.0d0))
    call BlockMat_set_element(mat, 1, 1, 1, 2, (1.0d0, 0.1d0))
    call BlockMat_set_element(mat, 1, 1, 2, 1, (1.0d0, 0.2d0))
    call BlockMat_set_element(mat, 1, 1, 2, 2, (1.0d0, 0.3d0))
    call BlockMat_set_element(mat, 2, 2, 1, 1, (2.0d0, 0.0d0))
    
    call BlockMat_set_element(mat, 3, 3, 1, 2, (3.0d0, 0.1d0))
    call BlockMat_set_element(mat, 3, 3, 2, 1, (3.0d0, 0.2d0))
    call BlockMat_set_element(mat, 3, 3, 2, 2, (3.0d0, 0.3d0))

    call BlockMat_set_element(mat, 1, 2, 1, 1, (5.0d0, 0.0d0))
    call BlockMat_set_element(mat, 1, 2, 1, 2, (5.0d0, 0.1d0))

    call BlockMat_block(mat, 1, 1, m11)
    call Expect_Eq("m11_11", (1.0d0, 0.0d0), m11(1, 1))
    call Expect_Eq("m11_12", (1.0d0, 0.1d0), m11(1, 2))
    call Expect_Eq("m11_21", (1.0d0, 0.2d0), m11(2, 1))
    call Expect_Eq("m11_22", (1.0d0, 0.3d0), m11(2, 2))

    call BlockMat_block(mat, 1, 2, m12)
    call Expect_Eq("m12_11", (5.0d0, 0.0d0), m12(1, 1))
    call Expect_Eq("m12_21", (5.0d0, 0.1d0), m12(2, 1))
    
  end subroutine test_BlockMat_get_as_array
  subroutine test_AoInts()
    type(AoInts) :: ao_ints
    integer, parameter :: ifile = 13
    open(unit=ifile, file='AOINTS', status='old', form='unformatted')
    call AoInts_new_read(ao_ints, ifile)
    call AoInts_show(ao_ints)
    call AoInts_delete(ao_ints)
    close(ifile)
  end subroutine test_AoInts
  subroutine test_DIPINTS()
    type(AoInts) :: ao_ints
    integer, parameter :: ifile = 13
    open(unit=ifile, file='DIPINTS', status='old', form='unformatted')
    call AoInts_new_read(ao_ints, ifile)
    call AoInts_show(ao_ints)
    call AoInts_delete(ao_ints)
    close(ifile)
  end subroutine test_DIPINTS
 
end program main

