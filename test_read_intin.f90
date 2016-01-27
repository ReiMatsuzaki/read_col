program main
  use Mod_IntIn
  implicit none
  integer, parameter :: ifile = 3
  type(IntIn) int_in
  integer is, ic, isf, icon, if, igto, igcs, ics, ict, ir, idx_gto
  integer iaords, ist
  integer lmp1
  integer, allocatable :: idx_ist(:)
  complex*16, allocatable :: vals_ist(:,:)

  call IntIn_new_read(int_in)

  allocate(idx_ist(int_in % nst))
  idx_ist(:) = 0
  allocate(vals_ist(int_in % nst, 100))
  
  isf = 0
  !  is : symmetry distinct atom set
  do is = 1, int_in % ns
     !  if : functions
     do if = 1, int_in % nf(is)
        isf = isf + 1
        lmp1 = int_in % lmnp1(isf)
        igto = int_in % mcons(isf)   ! index of GTO
        igcs = int_in % mgcs(isf)    ! index of GCS
        
        iaords= int_in % maords(igcs)
        !  ir,ics : symmetry index
        do ir = 1, int_in % nr(iaords)  
           ics = ir
           ist = int_in % la(ir, iaords)
           write(*, '(I2, "  ")', advance='no') igto
           write(*, '(A4)', advance='no') int_in % ityp(ist)
           
           ict = 0
           do ic = 1, int_in % nc(is)
              do idx_gto = 1, num_gto(lmp1)
                 ict = ict + 1
                 write(*, '(f4.1)', advance='no') int_in % igc(ict, ics, igcs)
                 write(*, '("[",3I1,"]@",I1,A3)', advance='no') &
                      nx(lmp1, idx_gto), ny(lmp1, idx_gto), nz(lmp1, idx_gto), &
                      ic, int_in % mtype(is)
              end do
              
              do icon = 1, int_in % ncon(isf)
              end do
              
           end do
           write(*, *) ""
           
        end do
     end do
  end do

end program main
!subroutine practice_proj_wave(int_in)
!  type(int_in), intent(in) :: int_in
!end subroutine practice_proj_wave
