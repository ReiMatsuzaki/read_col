! compute projection of N-1 electron wave function to
! each N electron determinant to obtain one particle orbital

module Mod_ProjOne
  use Mod_CI
  use Mod_AoInts
  use Mod_MOCoef
  implicit none  
  character(32)      :: in_path
  type(AoInts), save :: ao
  type(MOCoef), save :: mo
  type(CI), save     :: ci_wf

contains
  ! ==== general ====
  subroutine run()
    
    if(iargc() .ne. 1) then
       call print_usage
    end if
    
    call getarg(1, in_path)
    write(*, *) "in_path: ", in_path
    
    call new_read()

    if(ci_wf % n_ele == 2) then
       call two_e_case()
    else
       write(*, *) "only two electron systems can be treated"
       stop
    end if
        
  end subroutine run
  subroutine new_read()
    integer, parameter :: ifile_in = 13
    integer, parameter :: ifile_aoints = 14
    integer, parameter :: ifile_mocoef = 15
    integer, parameter :: ifile_csf = 16
    integer, parameter :: ifile_ciin = 17
    character(32)      :: aoints_path  = "no_file"
    character(32)      :: mocoef_path = "no_file"
    character(32)      :: csf_path  = "no_file"
    character(32)      :: ciin_path = "no_file"    
    integer n_ist(8)

    namelist/proj_one/aoints_path, mocoef_path, csf_path, ciin_path

    open(unit = ifile_in, file = in_path, status='old')
    read(unit = ifile_in, nml = proj_one)
    close(unit = ifile_in)

    write(*, *) "aoints_path: ", aoints_path
    write(*, *) "mocoef_path: ", mocoef_path
    write(*, *) "ciin_path: ", ciin_path
    write(*, *) "csf_path: ", csf_path

    open(unit = ifile_aoints, file = aoints_path, status = 'old', form='unformatted')
    call AoInts_new_read(ao, ifile_aoints)
    close(unit = ifile_aoints)

    n_ist(:) = ao % nso(1 : ao % nst)
    call MOCoef_new(mo, n_ist)
!    open(unit = ifile_mocoef, file = mocoef_path, status = 'old')
    call MOCoef_set_read(mo)
 !   close(unit = ifile_aoints)    

    open(unit = ifile_ciin, file = ciin_path, status = 'old')
    call CI_new_read_ciin(ci_wf, ifile_ciin)
    call CI_new_read_civec(ci_wf)
    close(unit = ifile_ciin)
  
    open(unit = ifile_csf, file = csf_path, status = 'old')
    call CI_new_read_csf(ci_wf, ifile_csf)
    close(unit = ifile_csf)
        
  end subroutine new_read
  subroutine print_usage()
    write(*, *) "proj_one INPUT"
    write(*, *) "------------------"
    write(*, *) "&proj_one"
    write(*, *) '"csf_path="CSF"'
    write(*, *) '"ciin_path="ci2.in"""'
    write(*, *) '"num_e=2' 
    write(*, *) "/"
    write(*, *) "------------------"
  end subroutine print_usage
  ! ==== two electron case ====
  subroutine two_e_case()
    use Mod_BlockVec
    complex*16, allocatable :: res(:, :, :)
    integer, parameter :: ifile = 13
    integer :: spin0 = 1 ! spin of final state ion
    integer :: isym0 = 1 ! symmetry of final state ion
    complex*16 :: ene0, x
    complex*16, allocatable :: cs0(:)
    character(32) :: symvec_path
    integer :: i_csf, i, i1, i2, mo1, mo2, mo1_sym, mo2_isym, ig1, ig2
    type(CSF) :: csf    
    type(BlockVec) :: block_vec, block_vec_mo

    namelist/two_e/spin0, isym0, symvec_path
    open(unit = ifile, file = in_path, status='old')
    read(unit = ifile, nml = two_e)
    close(unit = ifile)

    write(*, *)
    write(*, *) "Two electron case"
    write(*, *) "spin0: ", spin0
    write(*, *) "isym0: ", isym0
    write(*, *) "symvec_path: ", symvec_path

    write(*, *) "read BlockVec"
    open(unit = ifile, file = symvec_path)
    call BlockVec_new_read(block_vec, ifile)
    close(unit = ifile)

    write(*, *) "MO transform"
    call trans_symvec(block_vec, block_vec_mo)

    write(*, *) "diagonalize core hamiltonian"
    call diag_one_e(isym0, ene0, cs0)
    write(*, *) "E0: ", ene0
    write(*, *) "C0: "
    do i = 1, size(cs0)
       write(*, *) cs0(i)
    end do
    
    allocate(res(ci_wf % ncidim, 2, 2))
    
    do i_csf = 1, ci_wf % ncidim
       csf = ci_wf % csf(i_csf)
       do i_sd = 1, csf % n_sd
          csf % coef(i_sd)
          do i1 = 1, 2
             i2 = cond(i1 == 1, 2, 1)
             mo1 = csf % mo(i_sd, i1)
             mo1_sym = csf % mo_sym(i_sd, i1)
             mo2 = csf % mo(i_sd, i2)
             mo2_sym = csf % mo_sym(i_sd, i2)

             x = (0.0d0, 0.0d0)
             if(isym0 == mo2_sym) then
                do ig1 = 1, 2
                   do ig2 = 1, 2
                      x = x + BlockVec_val(block_vec_mo, mo1_sym, mo1) * &
                           cs0(ig1) * &
                           BlockMat_val(oa % s_mat, mo2_sym, mo2_sym, ig1, ig2) * &
                           BlockMat_val(mo % mo_coef, mo2_sym, mo2_sym, ig2, mo2)
                   end do
                end do
             end if
          end do
       end do
    end do
    
  end subroutine two_e_case
  subroutine diag_one_e(isym0, ene0, cs0)
    use Mod_LinearAlgebra
    integer, intent(in) :: isym0
    complex*16, intent(out) :: ene0
    complex*16, allocatable, intent(out) :: cs0(:)
    complex*16, allocatable :: s(:, :), h(:, :), v(:, :)
    complex*16, allocatable :: eigs(:), vecs(:, :)
    integer :: num_ij(2), n

    ! ==== size ====
    call BlockMat_check_block(ao % s_mat, isym0, isym0, "s_mat")
    call BlockMat_check_block(ao % t_mat, isym0, isym0, "t_mat")
    call BlockMat_check_block(ao % v_mat, isym0, isym0, "v_mat")
    num_ij(:) = BlockMat_block_size(ao % s_mat, isym0, isym0)
    if(num_ij(1) .ne. num_ij(2)) then
       write(*, *) "Error"
       write(*, *) "matrix is not squared"
       stop
    end if
    n = num_ij(1)

    ! ==== allocation ====
    allocate(s(n, n)); allocate(h(n, n)); allocate(v(n, n));
    allocate(eigs(n)); allocate(vecs(n, n))
    allocate(cs0(n))

    ! ==== set matrix ====
    call BlockMat_block(ao % s_mat, isym0, isym0, s)
    call BlockMat_block(ao % t_mat, isym0, isym0, h)
    call BlockMat_block(ao % v_mat, isym0, isym0, v)
    h(:, :) = h(:, :) + v(:, :)

    ! ==== solve ====
    call eigen_sym(h, s, eigs, vecs)
    ene0 = eigs(1)
    cs0 = vecs(:, 1)

    ! ==== free ====
    deallocate(s);    deallocate(h); deallocate(v);
    deallocate(eigs); deallocate(vecs)
    
  end subroutine diag_one_e
  subroutine trans_symvec(vec_in, vec_out)
    ! transform gto basis symmetry separated vector to MO basis
    type(BlockVec), intent(in) :: vec_in
    type(BlockVec), intent(out) :: vec_out

    vec_out = block_matmul(mo % mo_coef, vec_in, 't')
    
  end subroutine trans_symvec
  subroutine delta_lm(L, M, r0, res)
    integer, intent(in) :: L, M
    real*8,  intent(in) :: r0
    
    complex*16, allocatable, intent(out) :: res
    
  end subroutine delta_lm
end module Mod_ProjOne


