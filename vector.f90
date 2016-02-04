module Mod_VectorI
  type VectorI
     integer, allocatable :: data(:)
     integer :: num
  end type VectorI
contains
  function VectorI_new(num)
    type(VectorI) VectorI_new
    integer, intent(in), optional :: num
    if(present(num)) then
       allocate(VectorI_new % data(num))
    else
       allocate(VectorI_new % data(1))
    end if
    VectorI_new % num = 0
  end function VectorI_new
  subroutine VectorI_append(this, x)
    type(VectorI), intent(inout) :: this
    integer num_new
    integer, intent(in) :: x
    integer, allocatable :: tmp(:)
    if(size(this % data) == this % num) then
       num_new = this % num *2
       allocate(tmp(num_new))
       tmp(1:this % num) = this % data(:)
       deallocate(this % data)
       allocate(this % data(num_new))
       this % data(:) = tmp(:)
    end if

    this % data(this % num + 1) = x
    this % num = this % num + 1
    
  end subroutine VectorI_append
  function VectorI_at(this, i)
    type(VectorI) :: this 
    integer, intent(in) :: i
    integer :: VectorI_at
    VectorI_at = this % data(i)
  end function VectorI_at
  function VectorI_size(this)
    type(VectorI) :: this 
    integer VectorI_size
    VectorI_size = this % num
  end function VectorI_size
end module Mod_VectorI
module Mod_VectorZ
  type VectorZ
     complex*16, allocatable :: data(:)
     integer :: num
  end type VectorZ
contains
  function VectorZ_new(num)
    type(VectorZ) VectorZ_new
    integer, intent(in), optional :: num
    if(present(num)) then
       allocate(VectorZ_new % data(num))
    else
       allocate(VectorZ_new % data(1))
    end if
    VectorZ_new % num = 0
  end function VectorZ_new
  subroutine VectorZ_append(this, x)
    type(VectorZ), intent(inout) :: this
    integer num_new
    complex*16, intent(in) :: x
    complex*16, allocatable :: tmp(:)
    if(size(this % data) == this % num) then
       num_new = this % num *2
       allocate(tmp(num_new))
       tmp(1:this % num) = this % data(:)
       deallocate(this % data)
       allocate(this % data(num_new))
       this % data(:) = tmp(:)
    end if

    this % data(this % num + 1) = x
    this % num = this % num + 1
    
  end subroutine VectorZ_append
  function VectorZ_at(this, i)
    type(VectorZ) :: this 
    integer, intent(in) :: i
    complex*16 :: VectorZ_at
    VectorZ_at = this % data(i)
  end function VectorZ_at
  function VectorZ_size(this)
    type(VectorZ) :: this 
    integer VectorZ_size
    VectorZ_size = this % num
  end function VectorZ_size
end module Mod_VectorZ
module Mod_VectorD
  type VectorD
     real*8, allocatable :: data(:)
     integer :: num
  end type VectorD
contains
  function VectorD_new(num)
    type(VectorD) VectorD_new
    integer, intent(in), optional :: num
    if(present(num)) then
       allocate(VectorD_new % data(num))
    else
       allocate(VectorD_new % data(1))
    end if
    VectorD_new % num = 0
  end function VectorD_new
  subroutine VectorD_append(this, x)
    type(VectorD), intent(inout) :: this
    integer num_new
    real*8, intent(in) :: x
    real*8, allocatable :: tmp(:)
    if(size(this % data) == this % num) then
       num_new = this % num *2
       allocate(tmp(num_new))
       tmp(1:this % num) = this % data(:)
       deallocate(this % data)
       allocate(this % data(num_new))
       this % data(:) = tmp(:)
    end if

    this % data(this % num + 1) = x
    this % num = this % num + 1
    
  end subroutine VectorD_append
  function VectorD_at(this, i)
    type(VectorD) :: this 
    integer, intent(in) :: i
    real*8 :: VectorD_at
    VectorD_at = this % data(i)
  end function VectorD_at
  function VectorD_size(this)
    type(VectorD) :: this 
    integer VectorD_size
    VectorD_size = this % num
  end function VectorD_size
end module Mod_VectorD
