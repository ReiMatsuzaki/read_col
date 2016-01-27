module Mod_IntIn
  implicit none
  integer :: nx(5,10)
  integer :: ny(5,10)
  integer :: nz(5,10)
  integer :: num_gto(5) = (/1, 3, 6, 10, 0/)
  !                     ! lmp1=1 : S : 1
  !                     ! lmp1=2 : P : 3
  !                     ! lmp1=3 : D : 6
  !                     ! lmp1=4 : F : 10
  
  type IntIn
     character*100 blabel
     real*8 cscale, thetad ! used for complex scaling
     integer notwo, negl, irstart, itlim, nolim
     integer ngen  ! # of automaticcaly generated 
     integer ns    ! # of symmetry distinct atoms
     integer naords ! # of AO ReDuction Sets
     integer ncons  ! # of CONtracted GTOs
     integer ngcs   ! # of GTO coefficient set
     integer itol, icut
     integer nst ! # of symmetry
     !     integer ndpt ! # of 
     integer nd(10) ! (nd(ist), ist=1,nst) is digeneracy of symmetry(ist )
     character*4 ityp(10) ! (ityp(ist),ist=1,nst) is name of symmetry (ist)
     integer idp(8, 8, 8) ! point group product table. (1=>exist product, 0=>none)

     ! ==== AO ReDuction Sets ====
     integer nr(10)       ! (nr(iaords), iaords<-naords) : # of symmetry in iaords
     integer la(10, 10)   ! (la(IR, IAORDS),IR=1,..,) = symmetry sets for IAORDS:

     ! ==== GTO Contractions ====
     ! IGCS <- NGCS
     integer maords(10)      ! maords(IGCS) : IAORDS
     integer nct(10)         ! nct(IGCS) : # of GTOs
     integer ncs(10)         ! ncs(ICS)  : # of symmetry (almost same as nr)
     real*8  igc(10, 10, 10) ! igc(ICT, ICS, IGCS) : square coefficient
     !                       ! ICT<-nct(IGCS) : index of GTO
     !        (ex(f). 1=>XXX,2=>YYY,3=>ZZZ,4=>XXY,5=>XXZ,6=>YYZ,... )
     !                       ! ICS<-ncs(IGCS) : index of symmetry
     !                       ! IGCS : index of GTO Reduction Sets

     ! ==== GTOs ====
     integer lmnp1(108)      ! (lmnp1(ICONS), ICONS<-NCONS) : L+1
     character*4 istar(108)  ! (istar(ICONS), ICONS<-NCONS) : not used
     integer ncon(100)       ! (ncon(ICONS),  ICONS<-NCONS) : number of contraction 
     complex*16 zet(10, 100) ! (zet(icon, icons)) orbital exponent
     complex*16 eta(10, 100) ! coefficient

     ! ==== Atoms ====
     character*4 mtype(10)   ! (mtype(IS), IS<-NS) : name of atom
     integer nf(10)          ! (nf(IS), IS<-NS) : # of basis function
     integer nc(10)          ! (nc(IS), IS<-NS) : # of atoms
     real*8 chg(10)          ! (chg(IS), IS<-NS) : charge
     real*8 x(10, 8)         ! (x(IC, IS), IS<-NS, IC<-NC(IS))
     real*8 y(10, 8)
     real*8 z(10, 8)
     integer ica(12, 10, 49) ! ICA(IC,IS,IG)

     ! ==== Symmetry Adapted Basis ====
     integer mcons(108)      ! (MCONS(ISF)) index of contracted GTO
     integer mgcs(108)       ! (MGCS(ISF)) index of GCS
     
  end type IntIn
  type SymGtos
     ! atoms
     integer :: ns
     integer:: nc_is(10)
     real*8 :: x_ic_is(1, 1)
     real*8 :: y_ic_is(1, 1)
     real*8 :: z_ic_is(1, 1)

     ! Symmetry
     integer     :: nst          ! number of symmetry
     character*4 :: ityp_ist(10) ! name of symmetry
     
     ! Symmetry adapted basis
     integer :: num_sab_ist(10) ! number of Symmetry Adapted Basis
     
  end type SymGtos
contains
  subroutine IntIn_init()

    ! S
    nx(1, :) = 0
    ny(1, :) = 0
    nz(1, :) = 0

    ! P
    nx(2, :) = (/1, 0, 0,   1, 0, 0, 1, 0, 0, 1/)
    ny(2, :) = (/0, 1, 0,   0, 1, 0, 0, 1, 0, 0/)
    nz(2, :) = (/0, 0, 1,   0, 0, 1, 0, 0, 1, 0/)

    ! D
    nx(3, :) = (/2, 0, 0, 1, 1, 0,   9, 9, 9, 9/)
    ny(3, :) = (/0, 2, 0, 1, 0, 1,   9, 9, 9, 9/)
    nz(3, :) = (/0, 0, 2, 0, 1, 1,   9, 9, 9, 9/)

    ! F
    nx(4, :) = (/3, 0, 0, 2, 2, 1, 0, 1, 0, 1/)
    ny(4, :) = (/0, 3, 0, 1, 0, 2, 2, 0, 1, 1/)
    nz(4, :) = (/0, 0, 3, 0, 1, 0, 1, 2, 2, 1/)
    
  end subroutine IntIn_init
  subroutine IntIn_new_read(this)
    type(IntIn) this
!    integer ifile

    integer ist, idpt, jst, kst, iaords, iru, ir, igcs, icsu, ictu, ics, ict, icons, iconu, icon, isf, is, ic, ig, ng, if, ndpt
    integer igc2(10), it
    call IntIn_init()
    
    !    read(*, "(1H1, 8A10)") this % blabel
    read(*, *) this % blabel    
    read(*, '((2F16.10,9I3))') this % cscale, this % thetad, this % notwo, this % negl
    read(*, '(26I3)') this % ngen, this % ns, this % naords, this % ncons, this % ngcs, this % itol, this % icut

    read(*, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)

    read(*, *) ndpt
    do idpt = 1, ndpt
       read(*, *) ist, jst, kst
       this % idp(ist, jst, kst) = 1
       this % idp(ist, kst, jst) = 1
       this % idp(jst, ist, kst) = 1
       this % idp(jst, kst, ist) = 1
       this % idp(ist, jst, kst) = 1
       this % idp(kst, ist, jst) = 1
       this % idp(kst, jst, ist) = 1
    end do

    do iaords = 1, this % naords
       read(*, '(10I3)') iru, (this % la(ir, iaords), ir = 1, iru)
       this % nr(iaords) = iru
    end do

    do igcs = 1, this % ngcs
       read(*, '(10I3)') icsu, ictu, this % maords(igcs)
       this % nct(igcs) = ictu         ! nct(IGCS) : # of GTOs
       this % ncs(igcs) = icsu         ! ncs(ICS)  : # of symmetry
       do ics = 1, icsu
          read(*, '(10I3)') (igc2(ict), ict=1, ictu)
          do ict = 1, ictu
!             write(*, '(10I3)') igcs, ics, ict, igc2(ict)
             it = igc2(ict)
             if(it .eq. 0) then
                this % igc(ict, ics, igcs) = 0.0d0
             else
                this % igc(ict, ics, igcs) = it/abs(it) * sqrt(abs(real(it)))
             end if
          end do
       end do
    end do

    do icons = 1, this % ncons
       read(*, '(2I3,A4)') iconu, this % lmnp1(icons), this %istar(icons)
       read(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, iconu)
       this % ncon(icons) = iconu
    end do

    isf = 0
    do is = 1, this % ns
       read(*, *) this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       read(*, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       ng = this % ngen + 1
       if(this % nc(is) > 1) then
          do ig = 2, ng
             read(*, *) (this % ica(ic, is, ig), ic = 1, this % nc(is))
          end do
       end if
       do if = 1, this % nf(is)
          isf = isf + 1
          read(*, *) this % mcons(isf), this % mgcs(isf)
       end do
    end do
    
  end subroutine IntIn_new_read
  subroutine _IntIn_find_harmonics(igc2, lmp1, L, M, n)
    integer, intent(in) :: igc2(:)
    integer, intent(in) :: lmp1
    integer, intent(out) :: L, M, n
    integer :: num, i
    
    num = size(igc2)

    if (lmp1 .eq. 1) then
       if(igc2(1:1) .eq. (/1/)) then
          L = 0; M = 0; n = 0
          return 
       end if
       
    else if(lmp1 .eq. 2) then
       if(igc2(1:3) .eq. (/1, 0, 0/))  then
          L = 1; M = 1; n = 1
          return
       else if(igc2(1:3) .eq. (/0, 1, 0/)) then
          L = 1; M =-1; n = 1
          return
       else if(igc2(1:3) .eq. (/0, 0, 1/)) then
          L = 1; M = 0; n = 1
          return           
       end if

       ! D orbital
       ! (XX, YY, ZZ, XY, XZ, YZ)
    else if(lmp1 .eq. 3) then
       if(igc2(1:6) .eq. (/-1, -1, 4, 0, 0, 0/)) then
          L = 2; M = 0; n = 2; return 
       else if(igc2(1:6) .eq. (/1, -1, 0, 0, 0, 0/)) then
          L = 2; M = 2; n = 2; return
       else if(igc2(1:6) .eq. (/1, 1, 1, 0, 0, 0/)) then
          L = 0; M = 0; n = 2; return
       else if(igc2(1:6) .eq. (/0, 0, 0, 1, 0, 0/)) then
          L = 2; M = -2; n = 2; return 
       else if(igc2(1:6) .eq. (/0, 0, 0, 0, 1, 0/)) then          
          L = 2; M = +1; n = 2; return 
       else if(igc2(1:6) .eq. (/0, 0, 0, 0, 0, 1/)) then
          L = 2; M = -1; n = 2; return 
       end if

       ! F orbital
       ! 
    else if(lmp1 .eq. 4) then
    end if

    write(*, *) "Failed to find harmonics :", igc2
    stop
    
    do i = 1, num
       
    end do
    
  end subroutine IntIn_find_harmonics
  subroutine IntIn_show(this)
    type(IntIn) this
    integer icons, isf, ic, icon, if, is, ist, ng, iaords, ics, ict, igcs, ir

    call IntIn_init()
    
    write(*, *) this % blabel
    write(*, *) this % cscale, this % thetad, this % notwo, this % negl
    write(*, *) this % ngen, this % ns, this % naords, this % ngcs, this % itol, this % icut
    write(*, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)

    write(*, *) ""
    write(*, *) "AO ReDuction Sets"
    do iaords = 1, this % naords
       write(*, '(10I3)') (this % la(ir, iaords), ir = 1, this % nr(iaords))
    end do

    write(*, *) ""
    write(*, *) "Gto Contraction Sets"
    do igcs = 1, this % ngcs
       write(*, '(10I3)') this % ncs(igcs), this % nct(igcs), this % maords(igcs)
       do ics = 1, this % ncs(igcs)
          write(*, '(10I3)') (this % igc(ict, ics, igcs), ict=1, this%nct(igcs))
       end do
    end do
    do icons = 1, this % ncons
       write(*, '(2I3,A4)') this % ncon(icons), this % lmnp1(icons), this %istar(icons)
       write(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, this % ncon(icons))
    end do
    
    isf = 0
    do is = 1, this % ns
       write(*, *) this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       write(*, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       ng = this % ngen + 1
       do if = 1, this % nf(is)
          isf = isf + 1
          write(*, *) this % mcons(isf), this % mgcs(isf)
       end do
    end do
  end subroutine IntIn_show
  subroutine project_SH_GTO(eta, zeta, nx, ny, nz, L, M, r, res)
    complex*16, intent(in) :: eta, zeta
    integer, intent(in)    :: nx, ny, nz, L, M
    real*8, intent(in)     :: r
    complex*16, intent(out) :: res
    complex*16 :: rad
    integer :: n
    real*8 :: pi = atan(1.0d0)*4.0d0
    integer, parameter :: index_upper = 10
    integer, parameter :: nu = 3
    integer, parameter :: lu = 3
    complex * 16 :: val(index_upper)
    integer :: index(0:nu, 0:nu, 0:nu,  0:lu, -lu:lu)
    integer :: ix, iy, iz, l, m, idx
    complex*16, parameter :: ii = (0.0d0, 1.0d0)

    n = nx + ny + nz
    rad = eta * (r ** n)* exp(-zeta * r * r)
    index(:, :, :, :, :) = 0
    val(:) = (0.0d0, 0.0d0)
    idx = 0

    ! S
    idx = idx + 1; index(0, 0, 0, 0, 0) = idx
    val(idx) = rad * sqrt(4.0d0 * pi)

    ! Px
    idx = idx + 1; index(1, 0, 0, 1, 1) = idx
    val(idx) = -rad * sqrt(4.0d0 * pi / 3.0d0) / sqrt(2.0d0)
    idx = idx + 1; index(1, 0, 0, 1, -1) = idx
    val(idx) = +rad * sqrt(4.0d0 * pi / 3.0d0) / sqrt(2.0d0)
    ! Py
    idx = idx + 1; index(0, 1, 0, 1, 1) = idx
    val(idx) = -rad * sqrt(4.0d0 * pi / 3.0d0) / (sqrt(2.0d0) * ii)
    idx = idx + 1; index(0, 1, 0, 1,-1) = idx
    val(idx) = -rad * sqrt(4.0d0 * pi / 3.0d0) / (sqrt(2.0d0) * ii)
    ! Pz
    idx = idx + 1; index(0, 1, 0, 1, 0) = idx
    val(idx) = rad * sqrt(4.0d0 * pi / 3.0d0)

    ! D x2
    idx = idx + 1; index(2, 0, 0, 2, 0) = idx
    val(idx) = -0.5d0 * sqrt(16.0d0 * pi / 5.0d0) / 3.0d0 * rad
    idx = idx + 1; index(2, 0, 0, 0, 0) = idx
    val(idx) = +0.5d0 * 2.0d0 / 3.0d0 * sqrt(16.0d0 * pi / 5.0d0) * rad
    idx = idx + 1; index(2, 0, 0, 2, 2) = idx
    val(idx) = +0.5d0 * 1.0d0/sqrt(2.0d0) * sqrt(16.0d0 * pi / 15.0d0) * rad
    idx = idx + 1; index(2, 0, 0, 2,-2) = idx
    val(idx) = +0.5d0 * 1.0d0/sqrt(2.0d0) * sqrt(16.0d0 * pi / 15.0d0) * rad        
    ! D y2
    idx = idx + 1; index(0, 2, 0, 2, 0) = idx
    val(idx) = -0.5d0 * sqrt(16.0d0 * pi / 5.0d0) / 3.0d0 * rad
    idx = idx + 1; index(0, 2, 0, 0, 0) = idx
    val(idx) = +0.5d0 * 2.0d0 / 3.0d0 * sqrt(16.0d0 * pi / 5.0d0) * rad
    idx = idx + 1; index(0, 2, 0, 2, 2) = idx
    val(idx) = -0.5d0 * 1.0d0/sqrt(2.0d0) * sqrt(16.0d0 * pi / 15.0d0) * rad
    idx = idx + 1; index(2, 0, 0, 2,-2) = idx
    val(idx) = -0.5d0 * 1.0d0/sqrt(2.0d0) * sqrt(16.0d0 * pi / 15.0d0) * rad    
    ! D z2
    idx = idx + 1; index(0, 0, 2, 2, 0) = idx
    val(idx) = rad * sqrt(16.0d0 * pi / 5.0d0) / 3.0d0
    idx = idx + 1; index(0, 0, 2, 0, 0) = idx
    val(idx) = rad * sqrt(16.0d0 * pi / 4.0d0) / 3.0d0
    ! Dxy
    idx = idx + 1; index(1, 1, 0, 2, 2) = idx
    val(idx) = +rad * sqrt(4.0d0 * pi / 15.0d0) / (sqrt(2.0d0) * ii)
    idx = idx + 1; index(1, 1, 0, 2,-2) = idx
    val(idx) = -rad * sqrt(4.0d0 * pi / 15.0d0) / (sqrt(2.0d0) * ii)
    
    ! n = 1
    if(n .eq. 1) then
       if(L .ne. 1) then
          res = (0.0d0, 0.0d0)
          return
       else
          if(nx .eq. 1 .and. abs(M) .eq. 1) then
             res = - M * rad * sqrt(4.0d0 * pi / 3.0) * sqrt(1.0d0/2.0d0)
          else if(ny .eq. 1 .and. abs(M) .eq. 1) then
             res = - rad * sqrt(4.0d0 * pi / 3.0) * sqrt(1.0d0/2.0d0) / (0.0d0, 1.0d0)
          else if(nz .eq. 1 .and. M .eq. 0) then
             res = rad * sqrt(4.0d0 * pi / 3.0d0)
          else
             return             
          end if
       end if
    end if

    ! n = 2
    if(n .eq. 2) then
       if(  nx .eq. 1 .and. ny .eq. 1 .and. L .eq. 2) then
          if(M .eq. 1) then
             res = rad * sqrt(4.0d0*pi/15.0d0) / (sqrt(2.0d0) * (0.0d0, 1.0d0))
          else if(M .eq. -1) then
             res = -rad * sqrt(4.0d0*pi/15.0d0) / (sqrt(2.0d0) * (0.0d0, 1.0d0))
          end if
       else if( ny .eq. 1 .and. nz .eq. 1 .and. L .eq. 2) then
          if(M .eq. 1 .or. M .eq. -1) then
             res = -rad * sqrt(4.0d0*pi/15.0d0) / (sqrt(2.0d0) * (0.0d0, 1.0d0))
          else 
             return 
          end if          
       else if( nz .eq. 1 .and. nx .eq. 1 .and. L .eq. 2) then
          if(M .eq. 1) then
             res = -rad * sqrt(4.0d0*pi/15.0d0) / sqrt(2.0d0)
             return 
          else if(M .eq. -1) then
             res = +rad * sqrt(4.0d0*pi/15.0d0) / sqrt(2.0d0)
             return 
          end if          
       end if

    end if
    
    if(L .eq. 0 .and. M .eq. 0) then
       if(nx .eq. 0 .and. ny .eq. 0 .and. nz .eq. 0) then

       else
          res = (0.0d0, 0.0d0)
          return
       end if
    end if

    ! P
    
    
  end subroutine project_SH_GTO
  subroutine find_sh()
end module Mod_IntIn
