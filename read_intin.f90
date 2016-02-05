

!module Mod_GTO
!  implicit none
!  type GTO
!     integer :: num
!     integer :: lmp1
!     complex*16, allocatable :: zeta(:)
!     complex*16, allocatable :: eta(:)
!  end type GTO
!contains
!  subroutine GTO_new(this, lmp1, num)
!    type(GTO), intent(out) :: this
!    this % num = num
!    this % lmp1 = lmp1
!    allocate(this % zeta(num))
!    allocate(this % eta(num))
!  end subroutine GTO_new
!  subroutine GTO_delete(this)
!    type(GTO), intent(inout) :: this
!    deallocate(this % zeta)
!    deallocate(this % eta)
!  end subroutine GTO_delete
!  subroutine GTO_at(this, r, res)
!    type(GTO), intent(in) :: this
!    real*8, intent(in)    :: r
!    complex*16, intent(out)   :: res
!
!    integer i
!    res = (0.0d0, 0.0d0)
!    do i =  1, this % num
!       res = res + this % eta(i) * r ** this % lmp1 * exp(-this % zeta(i) * r * r)
!    end do
!  end subroutine GTO_at
!end module Mod_GTO

module Mod_GTO  
  use Mod_Utils
  implicit none
contains
  subroutine gto_int(pn, a, res)
    integer, intent(in) :: pn
    complex*16, intent(in) :: a
    complex*16, intent(out) :: res
    integer :: n, nn, k
    
    ! S_0^oo exp(-ax^2) dx = sqrt(pi/4a)
    ! S_0^oo x^2n exp(-ax^2) dx = (2n-1)! / (2 (2a)^n) x sqrt(pi/a)
    ! S_(-oo)^(oo) x^n exp(-zr^2) dx = 2

    if(pn .eq. 0) then
       res = sqrt(pi / (4.0d0 *a))
    else if(mod(pn, 2) .eq. 1) then
       n = (pn - 1) / 2
       nn = 1
       do k = 1, n
          nn = nn * k
       end do
       res = (1.0d0 * nn) / (2.0d0 * a**(n+1))
    else if(mod(pn, 2) .eq. 0) then
       n = pn / 2

       nn = 1
       do k = 1, 2*n-1, 2
          nn = nn * k
       end do       
       res = 1.0d0 * nn / (2.0d0 * (2.0d0*a) ** n)  * sqrt(pi / a)
    end if
    
  end subroutine gto_int
end module Mod_Gto

module Mod_IntIn
  use Mod_GTO
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
     integer igc2(10, 10, 10)
     real*8  igc(10, 10, 10) ! igc(ICT, ICS, IGCS) : square coefficient
     !                       ! ICT<-nct(IGCS) : index of GTO
     !        (ex(f). 1=>XXX,2=>YYY,3=>ZZZ,4=>XXY,5=>XXZ,6=>YYZ,... )
     !                       ! ICS<-ncs(IGCS) : index of symmetry
     !                       ! IGCS : index of GTO Reduction Sets
     integer :: gto_l_isf_ics(108, 10) ! L in real Spherical Harmonics
     integer :: gto_m_isf_ics(108, 10) ! M in real Spherical Harmonics
!     real*8 ::  gto_c_isf_ics(108, 10) ! coefficient of real Spherical Harmonics
     !                         ! GTO in real Spherical Harmonics form becomes
     !                              eta r^n exp(-zeta r^2) Y_LM

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

     ! ==== access from Symmetry Adapted Basis ====
     integer :: nsab_ist(8)
     integer :: isf_ist_isab(8, 108)
     integer :: is_ist_isab(8, 108)
     integer :: ics_ist_isab(8, 108)
     
  end type IntIn
  type SymGtos
     ! ==== Gtos ====
!     type(GTO), allocatable  :: gto_isf(:)
     integer :: L
     integer :: M
     real*16 :: coef_ylm

     ! ==== Atoms ====
     

     
     ! ==== Symmetry ====
     ! number of Symmetry
     ! ist <- nst
     integer :: nst  
     character*4, allocatable :: ityp_ist(:) ! name of symmetry

     ! (isf_list(i),i=offset_ist(ist),offset_ist(ist+1))
     ! gives index list of isf in ist symmetry     
     integer, allocatable :: isf_list(:)

     ! (ict_list(i),i=offset_ist(ist),offset_ist(ist+1))
     ! gives index list of ict in ist symmetry
     integer, allocatable :: ict_list(:)

     ! used for isf_list and ict_list
     integer, allocatable :: offset_ist(:)
     
     integer :: ns
     integer:: nc_is(10)
     real*8 :: x_ic_is(1, 1)
     real*8 :: y_ic_is(1, 1)
     real*8 :: z_ic_is(1, 1)
     
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
  subroutine IntIn_new_read(this, ifile, showq)
    type(IntIn), intent(inout) :: this
    integer, intent(in)        :: ifile
    logical, intent(in)        :: showq

    integer ist, idpt, jst, kst, iaords, iru, ir, igcs, icsu, ictu, ics, ict, icons, iconu, icon, isf, is, ic, ig, ng, if, ndpt
    integer it
    call IntIn_init()
    
    !    read(*, "(1H1, 8A10)") this % blabel
    read(ifile, *) this % blabel
    if(showq) then
       write(*, *) this % blabel
    end if
    read(ifile, '((2F16.10,9I3))') this % cscale, this % thetad, this % notwo, this % negl
    if(showq) then
       write(*, '((2F16.10,9I3))') this % cscale, this % thetad, this % notwo, this % negl
    end if
    read(ifile, '(26I3)') this % ngen, this % ns, this % naords, this % ncons, this % ngcs, this % itol, this % icut
    if(showq) then
       write(*, '(26I3)') this % ngen, this % ns, this % naords, this % ncons, this % ngcs, this % itol, this % icut
    end if
    

    read(ifile, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)
    if(showq) then
       write(*, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)
    end if

    ! ==== Group product table====
    read(ifile, *) ndpt
    if(showq) then
       write(*, *) ndpt
    end if
    do idpt = 1, ndpt
       read(ifile, *) ist, jst, kst
       if(showq) then
          write(*, *) ist, jst, kst
       end if
       this % idp(ist, jst, kst) = 1
       this % idp(ist, kst, jst) = 1
       this % idp(jst, ist, kst) = 1
       this % idp(jst, kst, ist) = 1
       this % idp(ist, jst, kst) = 1
       this % idp(kst, ist, jst) = 1
       this % idp(kst, jst, ist) = 1
    end do

    ! ==== AO ReDuction Sets ====
    do iaords = 1, this % naords
       read(ifile, '(10I3)') iru, (this % la(ir, iaords), ir = 1, iru)
       if(showq) then
          write(*, '(10I3)') iru, (this % la(ir, iaords), ir = 1, iru)
       end if
       this % nr(iaords) = iru
    end do
    
    ! ==== GTO contractions ====
    do igcs = 1, this % ngcs
       read(ifile, '(10I3)') icsu, ictu, this % maords(igcs)
       if(showq) then
          write(*, '(10I3)') icsu, ictu, this % maords(igcs)
       end if
       this % nct(igcs) = ictu         ! nct(IGCS) : # of GTOs
       this % ncs(igcs) = icsu         ! ncs(ICS)  : # of symmetry
       do ics = 1, icsu
          read(ifile, '(10I3)') (this % igc2(ict, ics, igcs), ict=1, ictu)
          if(showq) then
             write(*, '(10I3)') (this % igc2(ict, ics, igcs), ict=1, ictu)
          end if
          do ict = 1, ictu
             it = this % igc2(ict, ics, igcs)
             if(it .eq. 0) then
                this % igc(ict, ics, igcs) = 0.0d0
             else
                this % igc(ict, ics, igcs) = it/abs(it) * sqrt(abs(real(it)))
             end if
          end do
       end do
    end do

    ! ==== GTO ====
    do icons = 1, this % ncons
       read(ifile, '(2I3,A4)') iconu, this % lmnp1(icons), this %istar(icons)
       if(showq) then
          write(*, '(2I3,A4)') iconu, this % lmnp1(icons), this %istar(icons)
       end if
       read(ifile, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, iconu)
       if(showq) then
          write(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, iconu)
       end if
       this % ncon(icons) = iconu
    end do

    ! ==== Atoms ====
    isf = 0
    do is = 1, this % ns
       read(ifile, '(A3,2I3,F3.0)') this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       if(showq) then
          write(*, '(A3,2I3,F3.0)') this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       end if
       read(ifile, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       if(showq) then
          write(*, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       end if
       ng = this % ngen + 1
       if(this % nc(is) > 1) then
          do ig = 2, ng
             read(ifile, *) (this % ica(ic, is, ig), ic = 1, this % nc(is))
             if(showq) then
                write(*, *) (this % ica(ic, is, ig), ic = 1, this % nc(is))
             end if
          end do
       end if
       do if = 1, this % nf(is)
          isf = isf + 1
          read(ifile, *) this % mcons(isf), this % mgcs(isf)
          if(showq) then
             write(*, *) this % mcons(isf), this % mgcs(isf)
          end if
       end do
    end do


    call IntIn_set_spherical_harmonics(this)
    call IntIn_set_symmetry_access(this)
    call IntIn_normalize(this)
    
  end subroutine IntIn_new_read
  subroutine IntIn_set_spherical_harmonics(this)

    type(IntIn), intent(inout) :: this
    integer isf, is, if, icons, igcs, ics
    
    isf = 0
    do is = 1, this % ns
       do if = 1, this % nf(is)
          isf = isf + 1
          icons = this % mcons(isf)
          igcs  = this % mgcs(isf)
          do ics = 1, this % ncs(igcs)
             
             call IntIn_find_harmonics( &
                  this % igc2(:, ics, igcs), &
                  this % lmnp1(this % mcons(isf)), &
                  this % gto_l_isf_ics(isf, ics), &
                  this % gto_m_isf_ics(isf, ics))
          end do
       end do
    end do    
  end subroutine IntIn_set_spherical_harmonics
  subroutine IntIn_set_symmetry_access(this)
    type(IntIn), intent(inout) :: this
    integer isf, is, if, icons, igcs, iaords, ics, ist
    this % nsab_ist(:) = 0
    isf = 0
    do is = 1, this % ns
       do if = 1, this % nf(is)
          isf = isf + 1
          icons = this % mcons(isf)
          igcs  = this % mgcs(isf)
          iaords= this % maords(igcs)
          do ics = 1, this % ncs(igcs)
             ist = this % la(ics, iaords)
             this % nsab_ist(ist) = this % nsab_ist(ist) + 1
             this % isf_ist_isab(ist, this % nsab_ist(ist)) = isf
             this % is_ist_isab(ist, this % nsab_ist(ist)) = is
             this % ics_ist_isab(ist, this % nsab_ist(ist)) = ics
          end do
       end do
    end do
  end subroutine IntIn_set_symmetry_access
  subroutine IntIn_normalize(this)
    type(IntIn), intent(inout) :: this
    integer :: icons, icon, jcon
    complex*16 :: cumsum, norm, int_val, zet, eta
    integer :: pn

    ! u(r) = r^(n+1) Y_LM(^r) \sum_i eta_i exp(-zet_i r^2)
    ! (u, u) = \sum_ij eta_i eta_j (r^(2n+2) exp(-(zet_i+zet_j) r^2))

    
    do icons = 1, this % ncons
       pn = this % lmnp1(icons)
       cumsum = (0.0d0, 0.0d0)
       do icon = 1, this % ncon(icons)
          do jcon = 1, this % ncon(icons)
             zet = this % zet(icon, icons) * this % zet(jcon, icons)
             eta = this % eta(icon, icons) * this % eta(jcon, icons)
             call gto_int(2 * pn, zet, int_val)
             cumsum = cumsum + eta * int_val
          end do
       end do
       norm = sqrt(cumsum)
       do icon = 1, this % ncon(icons)
          this % eta(icon, icons) = this % eta(icon, icons) / norm
       end do
    end do
  end subroutine IntIn_normalize
  subroutine IntIn_find_harmonics(igc2, lmp1, L, M)

    ! x^n1 y^n2 z^n3 = coef r^n Y_LM
    
    integer, intent(in) :: igc2(:)
    integer, intent(in) :: lmp1
    integer, intent(out) :: L, M
    integer :: num
    real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
    real*8 :: c
    
    num = size(igc2)

    if (lmp1 .eq. 1) then
       if(all(igc2(1:1) .eq. (/1/))) then
          L = 0; M = 0; c = sqrt(4.0d0 * pi); return
       end if
       
    else if(lmp1 .eq. 2) then
       if(all(igc2(1:3) .eq. (/1, 0, 0/)))  then
          L = 1; M = 1; c = sqrt(4.0d0*pi/3.0d0); return
          return
       else if(all(igc2(1:3) .eq. (/0, 1, 0/))) then
          L = 1; M =-1; c = sqrt(4.0d0*pi/3.0d0); return
          return
       else if(all(igc2(1:3) .eq. (/0, 0, 1/))) then
          L = 1; M = 0; c = sqrt(4.0d0*pi/3.0d0); return 
       end if

       ! D orbital
       ! (XX, YY, ZZ, XY, XZ, YZ)
    else if(lmp1 .eq. 3) then
       if(all(igc2(1:6) .eq. (/-1, -1, 4, 0, 0, 0/))) then
          L = 2; M = 0; c = sqrt(16.0d0*pi/5.0d0); return 
       else if(all(igc2(1:6) .eq. (/1, -1, 0, 0, 0, 0/))) then
          L = 2; M = 2; c = sqrt(16.0d0*pi/15.0d0); return
       else if(all(igc2(1:6) .eq. (/1, 1, 1, 0, 0, 0/))) then
          L = 0; M = 0; c = sqrt(4.0d0*pi); return
       else if(all(igc2(1:6) .eq. (/0, 0, 0, 1, 0, 0/))) then
          L = 2; M = -2; c = sqrt(4.0d0*pi/15.0d0); return 
       else if(all(igc2(1:6) .eq. (/0, 0, 0, 0, 1, 0/))) then          
          L = 2; M = +1; c = sqrt(4.0d0*pi/15.0d0); return 
       else if(all(igc2(1:6) .eq. (/0, 0, 0, 0, 0, 1/))) then
          L = 2; M = -1; c = sqrt(4.0d0*pi/15.0d0); return 
       end if

       ! F orbital
       !                      ! (XXX,YYY,ZZZ,XXY,XXZ,YYX,YYZ,ZZX,ZZY,XYZ)
    else if(lmp1 .eq. 4) then
       ! 2z^3 -3xxz -3yyz = 5zzz - 3rrz
       if(all(igc2 .eq.      (/ 0,  0,  4,  0, -9,  0, -9,  0,  0,  0/))) then
          L = 3; M = 0; c = sqrt(16.0d0*pi/7.0d0); return 
       else if(all(igc2(1:10) .eq. (/-1,  0,  0,  0,  0, -1,  0, 16,  0,  0/))) then
          L = 3; M = 1; c = -sqrt(128.0d0*pi/21.0d0); return 
       else if(all(igc2(1:10) .eq. (/ 0, -1,  0, -1,  0,  0,  0,  0, 16,  0/))) then
          L = 3; M =-1; c = -sqrt(32.0d0*pi/21.0d0); return 
       else if(all(igc2(1:10) .eq. (/ 0,  0,  0,  0,  1,  0, -1,  0,  0,  0/))) then
          L = 3; M = 2; c = +sqrt(16.0d0*pi/105.0d0); return 
       else if(all(igc2(1:10) .eq. (/ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1/))) then
          L = 2; M = -2; c = -sqrt(16.0d0*pi/105.0d0); return 
       else if(all(igc2(1:10) .eq. (/ 0, -1,  0,  9,  0,  0,  0,  0,  0,  0/))) then
          L = 2; M = -3; c = -sqrt(32.0d0*pi/105.0d0); return 
       else if(all(igc2(1:10) .eq. (/-1,  0,  0,  0,  0,  9,  0,  0,  0,  0/))) then
          L = 2; M = +3; c = +sqrt(32.0d0*pi/105.0d0); return 
       end if
    end if

    write(*, *) "Failed to find harmonics :"
    write(*, *) "L+1", lmp1
    write(*, *) "igc2", igc2
    stop
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
          write(*, '(10f6.2)') (this % igc(ict, ics, igcs), ict=1, this%nct(igcs))
       end do
    end do
    do icons = 1, this % ncons
       write(*, '(2I3,A4)') this % ncon(icons), this % lmnp1(icons), this %istar(icons)
       !write(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, this % ncon(icons))
       write(*, *) (this % zet(icon, icons), this % eta(icon, icons), icon = 1, this % ncon(icons))
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
  subroutine IntIn_show_basis_symmetry(this)
    type(IntIn) :: this
    integer     :: ist, isab
    integer     :: isf, is, ics, ic, icon, ict, idx_gto, igcs, lmp1, icons

    do ist = 1, this % nst
       write(*, '(A4)') this % ityp(ist)
       do isab = 1, this % nsab_ist(ist)
          isf = this % isf_ist_isab(ist, isab)
          is  = this % is_ist_isab(ist, isab)
          ics = this % ics_ist_isab(ist, isab)
          icons = this % mcons(isf)
          igcs = this % mgcs(isf)
          lmp1 = this % lmnp1(icons)
          write(*, '(I2)', advance='no') this % mcons(isf)
          write(*, '(" Y(", I2, I2, ")" )', advance='no') &
               this % gto_l_isf_ics(isf, ics), &
               this % gto_m_isf_ics(isf, ics)


          ict = 0
          do ic = 1, this % nc(is)
             do idx_gto = 1, num_gto(lmp1)
                ict = ict + 1
                if(abs(this % igc(ict, ics, igcs)) > 0.001) then
                   write(*, '(f5.1)', advance='no') this % igc(ict, ics, igcs)
                   write(*, '("[",3I1,"]@",I1,A3)', advance='no') &
                        nx(lmp1, idx_gto), ny(lmp1, idx_gto), nz(lmp1, idx_gto), &
                        ic, this % mtype(is)
                end if
             end do
          end do
          
          write(*, *)
          do icon = 1, this % ncon(isf)
             write(*, '(4f10.5)') this % zet(icon, icons), this % eta(icon, icons)
          end do
          write(*, *) ""

       end do
    end do
    
  end subroutine IntIn_show_basis_symmetry
  subroutine IntIn_show_basis_symmetry2(this)
    type(IntIn) this
    integer is, ic, isf, icon, if, igto, igcs, ics, ict, ir, idx_gto
    integer iaords, ist
    integer lmp1
    
    isf = 0
    !  is : symmetry distinct atom set
    do is = 1, this % ns
       !  if : functions
       do if = 1, this % nf(is)
          isf = isf + 1
          lmp1 = this % lmnp1(isf)
          igto = this % mcons(isf)   ! index of GTO
          igcs = this % mgcs(isf)    ! index of GCS
          
          iaords= this % maords(igcs)
          !  ir,ics : symmetry index
          do ir = 1, this % nr(iaords)  
             ics = ir
             ist = this % la(ir, iaords)
             write(*, '(I2, "  ")', advance='no') igto
             write(*, '(A4)', advance='no') this % ityp(ist) 
             write(*, '("Y(", I2, I2, ")" )', advance='no') &
                  this % gto_l_isf_ics(isf, ics), &
                  this % gto_m_isf_ics(isf, ics)
             
             ict = 0
             do ic = 1, this % nc(is)
                do idx_gto = 1, num_gto(lmp1)
                   ict = ict + 1
                   if(abs(this % igc(ict, ics, igcs)) > 0.001) then
                      write(*, '(f5.1)', advance='no') this % igc(ict, ics, igcs)
                      write(*, '("[",3I1,"]@",I1,A3)', advance='no') &
                           nx(lmp1, idx_gto), ny(lmp1, idx_gto), nz(lmp1, idx_gto), &
                           ic, this % mtype(is)
                   end if
                end do
                
                do icon = 1, this % ncon(isf)
                end do
                
             end do
             write(*, *) ""
             
          end do
       end do
    end do
    
  end subroutine IntIn_show_basis_symmetry2
end module Mod_IntIn
