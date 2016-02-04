! Diagonalize Core hamiltonian and write its eigen values
! and eigen vectors.
!
! Input file
! ----------
! aoints : path for AOINTS (AOINTS)
! symmetry: calculation target symmetry index
! num_eigvec: # of eigen vectors which you want to obtain
!
! Output file
! (Input file infomation)
! eigen value 1
! eigen vector 1
! eigen value 2
! eigen vector 2
! ....


program main
  use Mod_AoInts
  use Mod_Utils
  use Mod_LinearAlgebra
  implicit none

  type(AoInts) :: ao
  integer, parameter :: in_file = 14
  integer, parameter :: aoints_file = 15
  integer :: isym, num_eig
  integer :: n, num(2), i, j
  character(100) :: in_path, aoints_path, dat_path
  complex*16, allocatable :: h_mat(:, :), s_mat(:, :), tmp_mat(:, :)
  complex*16, allocatable :: eigval(:), eigvec(:, :)

  namelist/solve_1e/isym,num_eig,aoints_path,dat_path

  if(iargc() .ne. 1) then
     call print_usage()
     stop
  end if

  ! ==== read input file ====
  call getarg(1, in_path)
  write(*, *) "in: ", in_path
  open(unit=in_file, file=in_path, status='old')
  read(unit=in_file, nml=solve_1e)

  write(*, *) "isym", isym
  write(*, *) "num_eig", num_eig
  write(*, *) "aoints_path", aoints_path
  write(*, *) "dat_path", dat_path

  close(unit = in_file)

  ! ==== read AOINTS file ====
  open(unit=aoints_file, file=aoints_path, status='old', form='unformatted')
  call AoInts_new_read(ao, aoints_file)
  close(unit = aoints_file)

  ! ==== memory allocation ====

  call SymBlockMat_check_block(ao % s_mat, isym, isym, "s_mat")
  call SymBlockMat_check_block(ao % t_mat, isym, isym, "t_mat")
  call SymBlockMat_check_block(ao % v_mat, isym, isym, "v_mat")
  
  num(:) = SymBlockMat_block_size(ao % s_mat, isym, isym)

  if(num(1) .ne. num(2)) then
     write(*, *) "matrix is not square"
     stop
  end if

  n = num(1)

  allocate(h_mat(n, n))
  allocate(s_mat(n, n))
  allocate(tmp_mat(n, n))
  allocate(eigval(n))
  allocate(eigvec(n, n))

  ! ==== set matrix and compute ====
  
  call SymBlockMat_block(ao % s_mat, isym, isym, s_mat)
  call SymBlockMat_block(ao % v_mat, isym, isym, h_mat)
  call SymBlockMat_block(ao % t_mat, isym, isym, tmp_mat)
  h_mat(:, :) = h_mat(:, :) + tmp_mat(:, :)

  call eigen_sym(h_mat, s_mat, eigval, eigvec)

  ! ==== write results ====
  write(*, *) isym
  write(*, *) n
  write(*, *) num_eig
  do i = 1, num_eig
     write(*, *) i, eigval(i)
     do j = 1, n
        write(*, *) eigvec(j, i)
     end do
  end do

  ! ==== Finalize ====
  call AoInts_delete(ao)
  deallocate(h_mat, s_mat, eigval, eigvec, tmp_mat)
  
contains
  subroutine print_usage()
    write(*, *) "usage: solve_1e INPUT"
    write(*, *) "   INPUT : name list input file for solve_1e"
    write(*, *) ""
     write(*, *) "Example of SOLVE_1E_IN:"
     write(*, *) "---------------------"
     write(*, *) "&solve_1e"
     write(*, *) "isym=1"
     write(*, *) "num_eig=1"
     write(*, *) 'aoints_path="AOINTS"'
     write(*, *) 'dat_path="solve_1e.dat"'
     write(*, *) "/"
     write(*, *) "---------------------"        
  end subroutine print_usage
end program main
