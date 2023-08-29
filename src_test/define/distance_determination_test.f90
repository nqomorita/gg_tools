!> 距離判定テストモジュール
module mod_gg_tools_distance_determination_test
  use mod_monolis_utils
  use mod_gg_tools
  implicit none

contains

  subroutine C3D8_shapefunc(n_base, local, func)
    implicit none
    integer(4), intent(in) :: n_base
    real(8), intent(in) :: local(3)
    real(8), intent(out) :: func(n_base)

    func(1) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0-local(3))
    func(2) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0-local(3))
    func(3) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0-local(3))
    func(4) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0-local(3))
    func(5) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0+local(3))
    func(6) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0+local(3))
    func(7) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0+local(3))
    func(8) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0+local(3))
  end subroutine C3D8_shapefunc

  subroutine C3D8_shapefunc_deriv(n_base, local, func)
    implicit none
    integer(kint), intent(in) :: n_base
    real(kdouble), intent(in) :: local(3)
    real(kdouble), intent(out) :: func(8,3)

    func(1,1) = -0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
    func(2,1) =  0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
    func(3,1) =  0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
    func(4,1) = -0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
    func(5,1) = -0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
    func(6,1) =  0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
    func(7,1) =  0.125d0*(1.0d0+local(2))*(1.0d0+local(3))
    func(8,1) = -0.125d0*(1.0d0+local(2))*(1.0d0+local(3))

    func(1,2) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
    func(2,2) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
    func(3,2) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
    func(4,2) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
    func(5,2) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(3))
    func(6,2) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
    func(7,2) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
    func(8,2) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(3))

    func(1,3) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(2,3) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(3,3) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(4,3) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
    func(5,3) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(6,3) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(7,3) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(8,3) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
  end subroutine C3D8_shapefunc_deriv

  subroutine gg_tools_distance_determination_test()
    implicit none
    integer(kint) :: n_base
    real(kdouble) :: local_pos(3), eps
    real(kdouble) :: deriv(8,3), ans(8,3)

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !>
        integer(4), intent(in) :: n_base
        !>
        real(8), intent(in) :: local_pos(3)
        !>
        real(8), intent(out) :: weight(n_base)
      end subroutine shape_func
    end interface

    procedure(shape_func), pointer ::  fptr => null()

    n_base = 8
    local_pos(1) = 0.0d0
    local_pos(2) = 0.0d0
    local_pos(3) = 0.0d0
    eps = 1.0d-3

    fptr => C3D8_shapefunc
    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)

    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    write(*,*)"ans"
    write(*,"(1p8e12.4)")deriv(:,1)
    write(*,"(1p8e12.4)")ans(:,1)
  end subroutine gg_tools_distance_determination_test

end module mod_gg_tools_distance_determination_test
