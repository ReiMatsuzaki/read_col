program main
  use Mod_AoInts
  use Mod_UnitTest  
  implicit none

  !  call test_SymBlockMat()
  !  call test_AoInts()
  call test_SymBlockMat()
  call test_SymBlockMat_get_as_array()
  
contains
  subroutine test_SymBlockMat()
    type(SymBlockMat) mat
    integer :: num_isym(4) = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    complex*16 val
    integer num_i, num_j
    logical have_val

    write(*, *) "test_SymBlockMat"
    call SymBlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    call IEq("iblock_ijsym(1,1)", 1, mat % iblock_ijsym(1, 1))
    call IEq("iblock_ijsym(2,2)", 2, mat % iblock_ijsym(2, 2))
    call IEq("iblock_ijsym(1,2)", 5, mat % iblock_ijsym(1, 2))
    call IEq("offset(1)", 0, mat % offset_iblock(1))
    call IEq("idx_1", 1, SymBlockMat_index(mat, 1, 1, 1, 1))

    call SymBlockMat_set(mat, 1, 1, 1, 1, (1.0d0, 0.0d0))
    call SymBlockMat_set(mat, 1, 1, 1, 2, (1.0d0, 0.1d0))
    call SymBlockMat_set(mat, 1, 1, 2, 1, (1.0d0, 0.2d0))
    call SymBlockMat_set(mat, 1, 1, 2, 2, (1.0d0, 0.3d0))

    call SymBlockMat_set(mat, 2, 2, 1, 1, (2.0d0, 0.0d0))
    
    call SymBlockMat_set(mat, 3, 3, 1, 2, (3.0d0, 0.1d0))
    call SymBlockMat_set(mat, 3, 3, 2, 1, (3.0d0, 0.2d0))
    call SymBlockMat_set(mat, 3, 3, 2, 2, (3.0d0, 0.3d0))

    call SymBlockMat_set(mat, 1, 2, 1, 1, (5.0d0, 0.0d0))
    call SymBlockMat_set(mat, 1, 2, 1, 2, (5.0d0, 0.1d0))

    ! call SymBlockMat_show(mat)
    call SymBlockMat_get(mat, 1, 1, 1, 1, val)
    call CEq("1,1,1,1", (1.0d0, 0.0d0), val)
    call SymBlockMat_get(mat, 1, 1, 1, 2, val)
    call CEq("1,1,1,2", (1.0d0, 0.1d0), val)
    call SymBlockMat_get(mat, 1, 2, 1, 2, val)
    call CEq("1,2,1,2", (5.0d0, 0.1d0), val)

    call SymBlockMat_get_block_size(mat, 1, 1, num_i, num_j ,have_val)
    call BTrue("block size(1,1)", have_val)
    call IEq("block size(1,1)", 2, num_i)
    call IEq("block size(1,1)", 2, num_j)

    call SymBlockMat_get_block_size(mat, 1, 3, num_i, num_j ,have_val)
    call BFalse("Block size(1,3)", have_val)

    call SymBlockMat_get_block_size(mat, 1, 2, num_i, num_j ,have_val)
    call BTrue("block size(1,2)", have_val)
    call IEq("block size(1,2)", 2, num_i)
    call IEq("block size(1,2)", 1, num_j)
    
    call SymBlockMat_delete(mat)
  end subroutine test_SymBlockMat
  subroutine test_SymBlockMat_get_as_array()
    type(SymBlockMat) mat
    integer :: num_isym(4) = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    complex*16 m11(2, 2), m12(2, 1)
    logical have_val

    write(*, *) "test_SymBlockMat"
    call SymBlockMat_new(mat, num_isym, isym_iblock, jsym_iblock)

    call SymBlockMat_set(mat, 1, 1, 1, 1, (1.0d0, 0.0d0))
    call SymBlockMat_set(mat, 1, 1, 1, 2, (1.0d0, 0.1d0))
    call SymBlockMat_set(mat, 1, 1, 2, 1, (1.0d0, 0.2d0))
    call SymBlockMat_set(mat, 1, 1, 2, 2, (1.0d0, 0.3d0))

    call SymBlockMat_set(mat, 2, 2, 1, 1, (2.0d0, 0.0d0))
    
    call SymBlockMat_set(mat, 3, 3, 1, 2, (3.0d0, 0.1d0))
    call SymBlockMat_set(mat, 3, 3, 2, 1, (3.0d0, 0.2d0))
    call SymBlockMat_set(mat, 3, 3, 2, 2, (3.0d0, 0.3d0))

    call SymBlockMat_set(mat, 1, 2, 1, 1, (5.0d0, 0.0d0))
    call SymBlockMat_set(mat, 1, 2, 1, 2, (5.0d0, 0.1d0))

    call SymBlockMat_get_as_array(mat, 1, 1, m11, have_val)
    call CEq("m11_11", (1.0d0, 0.0d0), m11(1, 1))
    call CEq("m11_12", (1.0d0, 0.1d0), m11(1, 2))
    call CEq("m11_21", (1.0d0, 0.2d0), m11(2, 1))
    call CEq("m11_22", (1.0d0, 0.3d0), m11(2, 2))

    call SymBlockMat_get_as_array(mat, 1, 2, m12, have_val)
    call CEq("m12_11", (5.0d0, 0.0d0), m12(1, 1))
    call CEq("m12_21", (5.0d0, 0.1d0), m12(2, 1))
    
  end subroutine test_SymBlockMat_get_as_array
  subroutine test_AoInts()
    type(AoInts) :: ao_ints
    integer, parameter :: ifile = 13
    open(unit=ifile, file='AOINTS', status='old', form='unformatted')
    call AoInts_new_read(ao_ints, ifile)
    call AoInts_show(ao_ints)
    call AoInts_delete(ao_ints)
    close(ifile)
  end subroutine test_AoInts
 
end program main

