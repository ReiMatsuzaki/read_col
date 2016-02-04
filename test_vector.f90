program main
  use Mod_UnitTest
  implicit none

  call run_test("VectorI", test_vectorI)
  call run_test("VectorZ", test_vectorZ)

contains
  subroutine test_vectorI()
    use Mod_VectorI
    type(VectorI) :: xs
    
    xs = VectorI_new(3)
    call VectorI_append(xs, 3)
    call VectorI_append(xs, 2)
    call VectorI_append(xs, 5)
    call VectorI_append(xs, 7)
    
    call expect_eq("xs(1)", 3, VectorI_at(xs, 1))
    call expect_eq("xs(2)", 2, VectorI_at(xs, 2))
    call expect_eq("size(xs)", 4, VectorI_size(xs))
  end subroutine test_vectorI
  subroutine test_vectorZ()
    use Mod_VectorZ
    type(VectorZ) :: xs

    xs = VectorZ_new()
    call VectorZ_append(xs, (1.1d0, 1.2d0))
    call VectorZ_append(xs, (1.1d0, 1.3d0))
    call VectorZ_append(xs, (1.1d0, 1.4d0))
    call VectorZ_append(xs, (1.1d0, 1.5d0))
    call expect_eq("xs(2)", (1.1d0, 1.3d0), VectorZ_at(xs, 2))
    call expect_eq("size(xs)", 4, VectorZ_size(xs))
  end subroutine test_vectorZ
end program main
