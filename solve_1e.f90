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
  implicit none

  type(AoInts) :: ao
  integer, parameter :: ifile = 13
  integer :: isym
  integer :: num_i, num_j, i, j
  complex*16, allocatable :: h_mat(:, :), s_mat(:, :)
  complex*16 :: val
  logical have_val
  complex*16, allocatable :: eigval(:), eigvec(:, :)

  isym = 1
  open(unit=ifile, file='AOINTS', status='old', form='unformatted')

  call AoInts_new_read(ao, ifile)

  call SymBlockMat_get_block_size(ao % s_mat, isym, isym, &
       num_i, num_j, have_val)

  if(.not. have_val) then
     write(*, *) "failed to find matrix for the symmetry"
  end if

  allocate(h_mat(num_i, num_j))
  allocate(s_mat(num_i, num_j))
  allocate(eigval(num_i))
  allocate(eigvec(num_i, num_j))

  do i = 1, num_i
     do j = 1, num_j
        call SymBlockMat_get(ao % s_mat, isym, isym, i, j, val)
        s_mat(i, j) = val
        call SymBlockMat_get(ao % t_mat, isym, isym, i, j, val)
        h_mat(i, j) = val
        call SymBlockMat_get(ao % v_mat, isym, isym, i, j, val)
        h_mat(i, j) = h_mat(i, j) + val
     end do
  end do

  call LA_GEEV_GEN(h_mat, s_mat, eigval, eigvec)

  ! ==== write results ====
  write(*, *) num_i
  write(*, *) isym
  do i = 1, num_i
     write(*, *) i, eigval(i)
     do j = 1, num_j
        write(*, *) eigvec(j, i)
     end do
  end do

  ! ==== Finalize ====
  call AoInts_delete(ao)
  deallocate(h_mat, s_mat, eigval, eigvec)
  close(unit=ifile)
  
end program main
