! Read CSF and CIVEC files produced by cCOLMBUS program and
! store CI coefficient, its eigene nergies and configurations.
module Mod_CI
  implicit none
  type CSF
     integer             :: n_sd    ! # of Slater Determinant
     integer             :: n_ele   ! # of electron (Dim of SD)
     real*8, allocatable :: coef(:) ! coef(i): coefficient of SD_i
     integer, allocatable:: mo_sym(:, :) ! mo_sym(i,j): symmetry of j th mo in SD_i
     integer, allocatable:: mo(:, :)! mo(i,j): j th mo in mo_sym(i,j) symmetry in SD_i
     integer, allocatable:: spin(:,:) ! 1=>alpha, -1=>beta
             
     !integer, allocatable:: mo(:, :)! mo(i, j): j th mo for SD_i
  end type CSF
  type CI
     integer   :: ncidim
     integer   :: nstate
     complex*16, allocatable :: eig(:)
!     complex*16              :: ene0
     complex*16, allocatable :: d_ene(:)
     complex*16, allocatable :: tdm_l(:)
     complex*16, allocatable :: tdm_v(:)
     complex*16, allocatable :: coef(:,:)
     integer, allocatable :: num_mo_sym(:)
     type(CSF), allocatable :: csf(:)
     integer, allocatable :: mo_offset(:)
  end type CI
  
contains
  subroutine CSF_new(this, n_sd, n_ele)
    type(CSF) this
    integer n_sd, n_ele
    this % n_sd = n_sd
    this % n_ele = n_ele
    allocate(this % coef(n_sd))
    allocate(this % mo_sym(n_sd, n_ele))
    allocate(this % mo(n_sd, n_ele))
    allocate(this % spin(n_sd, n_ele))
  end subroutine CSF_new
  subroutine CSF_delete(this)
    type(CSF) this
    deallocate(this % coef)
    deallocate(this % mo_sym)
    deallocate(this % mo)
  end subroutine CSF_delete
  subroutine CSF_show(this)
    type(CSF) this
    integer i, j
    write(*, *) ""
    write(*, *) "CSF Object"
    write(*, *) "n_sd", this % n_sd
    write(*, *) "n_ele", this % n_ele
    do i = 1, this % n_sd
       write(*, *) this % coef(i), (this % mo(i, j), this % mo_sym(i, j), j=1,this%n_ele)
    end do
  end subroutine CSF_show
  subroutine CI_new_read_ciin(this, ifile)
    
    type(CI) this
    integer ifile
    character*10 options
    integer i, num_mo, index, dum1, dum2
    integer, parameter :: num_sym = 8
    
    allocate(this % num_mo_sym(num_sym))
    allocate(this % mo_offset(num_sym))
    
    read(ifile, *) options
    do i = 1, num_sym
       read(ifile, *) num_mo, index, dum1, dum2
       this % num_mo_sym(i) = num_mo
    end do

    this % mo_offset(1) = 0
    do i = 2, num_sym
       this % mo_offset(i) = this % mo_offset(i-1) + this % num_mo_sym(i-1)
    end do
    
  end subroutine CI_new_read_ciin
  subroutine CI_new_read_civec(this)
    implicit none
    type(CI) this
    integer, parameter :: ifile = 12
    real*8 blabel(10) ! see ciprop.f
    integer name1(20) ! see ciprop.f
    integer ncidim, nstate, n

    open(unit=ifile, file='CIVEC', status='old', form='unformatted')

    read(ifile) ncidim, blabel, name1, nstate
    this % ncidim = ncidim
    this % nstate = nstate
    write(*, *) "ncidim: ", this % ncidim
    write(*, *) "nstate:   ", this % nstate

    allocate(this % eig(nstate))
    allocate(this % d_ene(nstate))
    allocate(this % tdm_l(nstate))
    allocate(this % tdm_v(nstate))
    allocate(this % coef(nstate,  ncidim))
    allocate(this % csf(ncidim))
        
    do n = 1, this % nstate
       read(ifile) this % eig(n), this % coef(n,:)
       !       write(*, *) this % eig(n)
    end do    
    close(unit=ifile)    
  end subroutine CI_new_read_civec
  subroutine CI_set_read_phoxsec(this, ifile)
    type(CI), intent(inout) :: this
    integer, intent(in)  :: ifile
    character(100) dum
    complex*16 :: ene, tdm_l, tdm_v
    integer i, ii
    
    do i = 1, 5
       read(ifile, *) dum
    end do

    read(ifile, '(I8,1P,6D12.4)') i, ene, tdm_l, tdm_v
!    this % ene0 = ene
    read(ifile, *) dum
    do i = 1, this % nstate
       read(ifile, '(I4,F12.7,5F12.8)') ii, ene, tdm_l, tdm_v
       this % d_ene(i) = ene
       this % tdm_l(i) = tdm_l
       this % tdm_v(i) = tdm_v
    end do
    
  end subroutine CI_set_read_phoxsec
  subroutine CI_set_read_csf(this, ifile)
    type(CI) this
    integer ifile
    !    integer :: num_mo_sym(:)    ! num_mo_sym(i) : # of MO for symmetry i
    character*8 label
    integer dum, n_ele, i, i_csf, i_sd, offset, i_sym, mo_i, spin_i
    integer, allocatable :: mo(:)
    real*8 coef

    read(ifile, *) label
    read(ifile, *) n_ele
    write(*, *) "n_ele:", n_ele
    allocate(mo(n_ele))
    i_csf = 1
    i_sd = 1
    do
       read(ifile, *, end=999) dum, coef, (mo(i),i=1,n_ele)
       if(dum .eq. 1 .or. dum .eq. 0) then
          this % csf(i_csf) % coef(i_sd) = coef
          do i = 1, n_ele
             mo_i = abs(mo(i))
             spin_i = mo_i / mo(i)
             do i_sym = 1, size(this % mo_offset)
                offset = this % mo_offset(i_sym)
                if(offset .lt. mo_i .and. offset .le. mo_i) then
                   this % csf(i_csf) % mo(i_sd, i)     = mo_i - offset
                   this % csf(i_csf) % mo_sym(i_sd, i) = i_sym
                   this % csf(i_csf) % spin(i_sd, i) = spin_i
                end if
             end do
          end do
       end if
       
       if(dum .eq. 1) then         ! exist next entry
          i_sd = i_sd + 1
       else if(dum .eq. 0) then
          i_csf = i_csf + 1
          i_sd = 1
       end if
       
    end do
999 continue
    deallocate(mo)
  end subroutine CI_set_read_csf
  subroutine CI_new_read_csf(this, ifile)
    type(CI) this
    integer ifile
    character*8 label
    integer dum, n_ele, i, i_csf, i_sd
    integer, allocatable :: mo(:)
    real*8 coef

    read(ifile, *) label
    read(ifile, *) n_ele
    write(*, *) "label: ", label
    write(*, *) "n_ele", n_ele
    allocate(mo(n_ele))
    i_csf = 0
    i_sd = 1
    do
       read(ifile, *, end=999) dum, coef, (mo(i),i=1,n_ele)
       if(dum .eq. 1) then ! exist next entry
          i_sd = i_sd + 1
       else if(dum .eq. 0) then
          i_csf = i_csf + 1
          call CSF_new(this % csf(i_csf), i_sd, n_ele)
          i_sd = 1
       else if(dum .eq. 3) then
       end if
    end do
999 continue    
    deallocate(mo)

    
    
  end subroutine CI_new_read_csf
  subroutine CI_delete(this)
    type(CI) this
    deallocate(this % coef)
    deallocate(this % eig)
  end subroutine CI_delete
end module Mod_CI


