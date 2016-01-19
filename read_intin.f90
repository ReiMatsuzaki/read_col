module Mod_IntIn
  implicit none
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
     character*4 ityp(10) ! ityp(ist) is name of symmetry ityp
     integer idp(8, 8, 8)
     integer la(10, 10)
     integer maords(10)
     integer igc(10, 10, 10)
     integer lmnp1(100)
     character*4 istar(100)
     complex*16 zet(10, 100)
     complex*16 eta(10, 100)
     character*4 mtype(10)
     integer nf(10)
     integer nc(10)
     real*8 chg(10)
     real*8 x(10, 8)
     real*8 y(10, 8)
     real*8 z(10, 8)
     integer ica(12, 10, 49)
     integer mcons(108)
     integer mgcs(108)
  end type IntIn
contains
  subroutine IntIn_new_read(this)
    type(IntIn) this
!    integer ifile

    integer ist, idpt, jst, kst, iaords, iru, ir, igcs, icsu, ictu, ics, ict, icons, iconu, icon, isf, is, ic, ig, ng, if, ndpt
    
    !    read(*, "(1H1, 8A10)") this % blabel
    read(*, *) this % blabel
    write(*, *) this % blabel
    
    read(*, '((2F16.10,9I3))') this % cscale, this % thetad, this % notwo, this % negl
    write(*, *) this % cscale, this % thetad, this % notwo, this % negl
    read(*, '(26I3)') this % ngen, this % ns, this % naords, this % ncons, this % ngcs, this % itol, this % icut
    write(*, *) this % ngen, this % ns, this % naords, this % ngcs, this % itol, this % icut
    read(*, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)
    
    
    write(*, '(I3, 12(I3,A3))') this % nst, (this % nd(ist), this % ityp(ist), ist=1, this % nst)

    read(*, *) ndpt
    write(*, *) ndpt
    do idpt = 1, ndpt
       read(*, *) ist, jst, kst
       this % idp(ist, jst, kst) = 1
       this % idp(ist, kst, jst) = 1
       this % idp(jst, ist, kst) = 1
       this % idp(jst, kst, ist) = 1
       this % idp(ist, jst, kst) = 1
       this % idp(kst, ist, jst) = 1
       this % idp(kst, jst, ist) = 1
       write(*, *) ist, jst, kst
    end do

    write(*, *)
    do iaords = 1, this % naords
       read(*, '(10I3)') iru, (this % la(ir, iaords), ir = 1, iru)
       write(*, '(10I3)') iru, (this % la(ir, iaords), ir = 1, iru)
    end do

    write(*, *)
    do igcs = 1, this % ngcs
       read(*, '(10I3)') icsu, ictu, this % maords(igcs)
       write(*, '(10I3)') icsu, ictu, this % maords(igcs)
       do ics = 1, icsu
          read(*, '(10I3)') (this % igc(ict, ics, igcs), ict=1, ictu)
          write(*, '(10I3)') (this % igc(ict, ics, igcs), ict=1, ictu)
       end do
    end do

    do icons = 1, this % ncons
       read(*, '(2I3,A4)') iconu, this % lmnp1(icons), this %istar(icons)
       read(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, iconu)
       write(*, '(2I3,A4)') iconu, this % lmnp1(icons), this %istar(icons)
       write(*, '(4D14.8)') (this % zet(icon, icons), this % eta(icon, icons), icon = 1, iconu)
    end do

    isf = 0
    do is = 1, this % ns
       read(*, *) this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       write(*, *) this % mtype(is), this % nf(is), this % nc(is), this % chg(is)
       read(*, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       write(*, *) (this % x(ic, is), this % y(ic, is), this % z(ic, is), ic = 1, this % nc(is))
       ng = this % ngen + 1
       if(this % nc(is) > 1) then
          do ig = 2, ng
             read(*, *) (this % ica(ic, is, ig), ic = 1, this % nc(is))
          end do
       end if
       do if = 1, this % nf(is)
          isf = isf + 1
          read(*, *) this % mcons(isf), this % mgcs(isf)
          write(*, *) this % mcons(isf), this % mgcs(isf)
       end do
    end do
    
  end subroutine IntIn_new_read
end module Mod_IntIn
