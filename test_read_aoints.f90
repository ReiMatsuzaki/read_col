program main
  use Mod_SymBlockMat
  !  use Mod_SymMat
  use Mod_UnitTest
  implicit none

  call test_SymBlockMat()
  
contains
  subroutine test_SymBlockMat()
    type(SymBlockMat) mat
    integer :: num_isym(4) = (/2, 1, 3, 4/)
    integer :: isym_iblock(5) = (/1, 2, 3, 4, 1/)
    integer :: jsym_iblock(5) = (/1, 2, 3, 4, 2/)
    complex*16 val
    integer num_i, num_j
    logical have_val

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

!  subroutine test_symmat()
!    type(SymMat) m1
!    real*8 eps
!    eps = 10.0**(-10.0)
!    
!    call SymMat_new(m1, 3)
!    
!    call SymMat_set(m1, 1, 1, (1.2d0, 0.0d0))
!    call SymMat_set(m1, 1, 2, (1.1d0, 0.0d0))
!    call SymMat_set(m1, 1, 3, (1.3d0, 0.0d0))
!    
!    call CNear("m1(1,2)", m1 % val(1, 3), (1.3d0, 0.0d0), eps)
!    call CNear("m1(2,1)", m1 % val(1, 2), m1 % val(2, 1), eps)
!      
!    call SymMat_delete(m1)
!  end subroutine test_symmat
!  subroutine test_sym_block_mat()
!    type(SymBlockMat) bmat
!    type(SymMat), pointer ::  mat1, mat2
!    logical non_zero
!
!    ! Initialize 
!    call SymBlockMat_new_diag(bmat, [2, 3, 4])
!
!    ! Check off diagonal is empty
!    call SymBlockMat_get_mat(bmat, 1, 2, mat1, non_zero)
!    call BFalse("(1,2) is empty ", non_zero)
!
!    ! Set (1, 1) block matrix
!    write(*, *) "A"
!    call SymBlockMat_get_mat(bmat, 1, 1, mat1, non_zero)
!    write(*, *) "B"
!    call BTrue("(1,1) is not empty ", non_zero)
!    call SymMat_set(mat1, 1, 1, (1.0d0, 1.1d0))
!    call SymMat_set(mat1, 1, 2, (1.0d0, 1.2d0))
!    call SymMat_set(mat1, 2, 2, (1.1d0, 1.2d0))
!
!    ! Set (2, 2) block matrix
!    call SymBlockMat_get_mat(bmat, 2, 2, mat1, non_zero)
!    call SymMat_set(mat1, 1, 1, (2.0d0, 1.0d0))
!    call SymMat_set(mat1, 1, 2, (2.0d0, 1.1d0))
!
!    ! Check (1, 1) block matrix
!    call SymBlockMat_get_mat(bmat, 1, 1, mat2, non_zero)
!    call CEq("bmat(1)(1,1)", (1.0d0, 1.1d0), mat2 % val(1, 1))
!    call CEq("bmat(1)(2,1)", (1.0d0, 1.2d0), mat2 % val(2, 1))
!
!    ! Check (2, 2) block matrix
!    call SymBlockMat_get_mat(bmat, 2, 2, mat2, non_zero)
!    call CEq("bmat(2,2)(1,1)", (2.0d0, 1.0d0)+0.1d0, mat2 % val(1, 1))
!    call CEq("bmat(2,2)(2,1)", (2.0d0, 1.1d0), mat2 % val(2, 1))
!
!    ! Show
!    call SymBlockMat_show(bmat)
!
!    ! Finalize
!    call SymBlockMat_delete(bmat)
!  end subroutine test_sym_block_mat
end program main

