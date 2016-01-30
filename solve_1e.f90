! Diagonalize Core hamiltonian and write its eigen values
! and eigen vectors.


program main
  use Mod_AoInts
  use f95_lapack, only : LA_GEES
  implicit none

  type(AoInts) :: ao
  integer, parameter :: ifile = 13
  integer :: isym
  integer :: num_i, num_j, i, j
  complex*16, allocatable :: h_mat(:, :), s_mat(:, :)
  complex*16 :: val
  logical have_val

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

  do i = 1, num_i
     do j = 1, num_j
        call SymBlockMat_get(ao % s_mat, isym, isym, i, j, val)
        s_mat(i, j) = val
        call SymBlockMat_get(ao % t_mat, isym, isym, i, j, val)
        h_mat(i, j) = -0.5d0 * val
        call SymBlockMat_get(ao % v_mat, isym, isym, i, j, val)
        h_mat(i, j) = h_mat(i, j) + val
     end do
  end do

  write(*, *) h_mat
  
  !  call AoInts_show(ao)
  call AoInts_delete(ao)
  
  close(unit=ifile)
  
  write(*, *) "hello"
  
end program main
