! This program calculate vector quantity using GTO.
! This need two input files. One is input file of molint
! and the other is for describing which quantity the
!user want to calculate.
! RM 2016/1/27

! name llist?

program main
  use Mod_SymVec
  use Mod_IntIn
  implicit none

  type(IntIn) int_in
  real*8 :: r0 = 7.7d0
  integer, parameter :: intin_file = 14
  integer, parameter :: in_file = 15
  integer, parameter :: dat_file = 13
  character(100) :: dat_path
  character(100) :: intin_path
  character(100) :: in_path

  namelist/gto_vec/r0, dat_path

  if(iargc() .ne. 2) then
     call print_usage()
     stop
  end if

  call getarg(1, intin_path)
  write(*, *) "int.in path: ", intin_path
  open(unit = intin_file, file=intin_path)
  call IntIn_new_read(int_in, intin_file, .false.)

  call getarg(2, in_path)
  write(*, *) "GTO_VEC_IN: ", in_path
  open(unit = in_file, file=in_path)
  read(unit = in_file, nml=gto_vec)
  
  open(unit = dat_file, file = dat_path)
  
  write(*, *) "r0: ", r0
  write(*, *) "dat_path: ", dat_path
  
  call calc_delta_ylm(int_in, r0, 1, 0, dat_file)
  
  close(unit = dat_file)
  close(unit = in_file)
  close(unit = intin_file)

contains
  subroutine print_usage()
     ! write(*, *) "gto_vec stopped. Number of argments must be 2"
     write(*, *) "usage: gto_vec MOLINT_IN GTO_VEC_IN"
     write(*, *) "    MOLINT_IN  : input file for molint."
     write(*, *) "    GTO_VEC_IN : name list input file for gto_vec."
     write(*, *) ""
     write(*, *) "Example of GTO_VEC_IN:"
     write(*, *) "---------------------"
     write(*, *) "&gto_vec"
     write(*, *) "r0=10.0"
     write(*, *) 'dat_path="gto_vec_delta.dat"'
     write(*, *) "/"
     write(*, *) "---------------------"    
  end subroutine print_usage
  subroutine calc_delta_ylm(int_in, r0, L, M, dat_file)
    type(IntIn), intent(in) :: int_in
    real*8, intent(in)      :: r0
    integer, intent(in)     :: L, M
    integer, intent(in)     :: dat_file

    type(SymVec) :: vec
    integer isf, is, ic, icon, igcs, ics
    integer icons, isab, ist
    integer L_y, M_y, lmp1
    complex*16 :: cumsum

    call SymVec_new(vec, int_in % nsab_ist(1:int_in % nst))

    do ist = 1, int_in % nst
       do isab = 1, int_in % nsab_ist(ist)
          isf = int_in % isf_ist_isab(ist, isab)
          is  = int_in % is_ist_isab(ist, isab)
          ics = int_in % ics_ist_isab(ist, isab)
          
          icons = int_in % mcons(isf)
          igcs = int_in % mgcs(isf)
          lmp1 = int_in % lmnp1(icons)
          L_y = int_in % gto_l_isf_ics(isf, ics)
          M_y = int_in % gto_m_isf_ics(isf, ics)

          if(L .eq. L_y .and. M .eq. M_y) then

             if(int_in % nc(is) .ne. 1) then
                write(*, *) "nc(is) != 1 is not supported"
                stop
             end if

             cumsum = (0.0d0, 0.0d0)
             do ic = 1, int_in % nc(is)
                do icon = 1, int_in % ncon(isf)
                   cumsum = cumsum + &
                        int_in % eta(icon, icons) * &
                        exp(-int_in % zet(icon, icons)*r0*r0)
                end do
             end do
             cumsum = cumsum * r0 ** lmp1 
             call SymVec_set(vec, ist, isab, cumsum)

          else
             call SymVec_set(vec, ist, isab, (0.0d0, 0.0d0))
          end if
          
       end do
    end do
    
    call SymVec_write(vec, dat_file)

  end subroutine calc_delta_ylm
end program main
