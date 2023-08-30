!> 距離判定テストモジュール
module mod_gg_tools_distance_determination_test
  use mod_monolis_utils
  use mod_gg_tools
  implicit none

contains

  subroutine gg_tools_distance_determination_test()
    implicit none

    call get_shapefunc_deriv_3d_test()
    call gg_tools_get_distance_3d_test()
    call gg_tools_get_distance_2d_test()
  end subroutine gg_tools_distance_determination_test

  subroutine gg_tools_get_distance_2d_test()
    implicit none
    integer(kint) :: n_base
    real(kdouble) :: coord(3,4)
    pointer :: shape_func
    real(kdouble) :: global_pos(3)
    real(kdouble) :: local_pos(2)
    real(kdouble) :: ths
    real(kdouble) :: ths_up
    real(kdouble) :: eps
    real(kdouble) :: ans(2)
    logical :: is_converge

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !> 要素を構成する形状関数の個数
        integer(4), intent(in) :: n_base
        !> 入力座標に最も近い局所座標
        real(8), intent(in) :: local_pos(2)
        !> 形状関数の重み
        real(8), intent(out) :: weight(n_base)
      end subroutine shape_func
    end interface

    procedure(shape_func), pointer ::  fptr => null()

    call monolis_std_global_log_string("gg_tools_get_distance_2d")

    n_base = 4
    eps = 1.0d-3
    ths = 1.0d-3
    ths_up = 1.0d+3

    fptr => C2D4_shapefunc

    coord(1,1) = 0.0d0; coord(2,1) = 0.0d0; coord(3,1) = 0.0d0
    coord(1,2) = 1.0d0; coord(2,2) = 0.0d0; coord(3,2) = 0.0d0
    coord(1,3) = 1.0d0; coord(2,3) = 1.0d0; coord(3,3) = 0.0d0
    coord(1,4) = 0.0d0; coord(2,4) = 1.0d0; coord(3,4) = 0.0d0

    !> case 1
    global_pos(1) = 0.0d0
    global_pos(2) = 0.0d0
    global_pos(3) = 0.0d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) =-1.0d0
    ans(2) =-1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 1", local_pos, ans)

    !> case 2
    global_pos(1) = 1.0d0
    global_pos(2) = 1.0d0
    global_pos(3) = 1.0d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 1.0d0
    ans(2) = 1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 2", local_pos, ans)

    !> case 3
    global_pos(1) = 0.5d0
    global_pos(2) = 0.5d0
    global_pos(3) = 0.5d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 0.0d0
    ans(2) = 0.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 3", local_pos, ans)

    !> pattern 2
    coord(1,1) = 0.0d0; coord(2,1) = 0.0d0; coord(3,1) = 0.0d0
    coord(1,2) = 1.0d0; coord(2,2) = 0.0d0; coord(3,2) = 0.0d0
    coord(1,3) = 1.0d0; coord(2,3) = 1.0d0; coord(3,3) = 0.0d0
    coord(1,4) = 0.0d0; coord(2,4) = 1.0d0; coord(3,4) = 0.0d0

    !> case 4
    global_pos(1) = 0.0d0
    global_pos(2) = 0.0d0
    global_pos(3) = 0.0d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) =-1.0d0
    ans(2) =-1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 4", local_pos, ans)

    !> case 5
    global_pos(1) = 1.0d0
    global_pos(2) = 1.0d0
    global_pos(3) = 1.0d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 1.0d0
    ans(2) = 1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 5", local_pos, ans)

    !> case 6
    global_pos(1) = 0.5d0
    global_pos(2) = 0.5d0
    global_pos(3) = 0.5d0

    call gg_tools_get_distance_2d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 0.0d0
    ans(2) = 0.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_2d_test 6", local_pos, ans)
  end subroutine gg_tools_get_distance_2d_test

  subroutine gg_tools_get_distance_3d_test()
    implicit none
    integer(kint) :: n_base
    real(kdouble) :: coord(3,8)
    pointer :: shape_func
    real(kdouble) :: global_pos(3)
    real(kdouble) :: local_pos(3)
    real(kdouble) :: ths
    real(kdouble) :: ths_up
    real(kdouble) :: eps
    real(kdouble) :: ans(3)
    logical :: is_converge

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !> 要素を構成する形状関数の個数
        integer(4), intent(in) :: n_base
        !> 入力座標に最も近い局所座標
        real(8), intent(in) :: local_pos(3)
        !> 形状関数の重み
        real(8), intent(out) :: weight(n_base)
      end subroutine shape_func
    end interface

    procedure(shape_func), pointer ::  fptr => null()

    call monolis_std_global_log_string("gg_tools_get_distance_3d")

    n_base = 8
    eps = 1.0d-3
    ths = 1.0d-3
    ths_up = 1.0d+3

    fptr => C3D8_shapefunc

    coord(1,1) = 0.0d0; coord(2,1) = 0.0d0; coord(3,1) = 0.0d0
    coord(1,2) = 1.0d0; coord(2,2) = 0.0d0; coord(3,2) = 0.0d0
    coord(1,3) = 1.0d0; coord(2,3) = 1.0d0; coord(3,3) = 0.0d0
    coord(1,4) = 0.0d0; coord(2,4) = 1.0d0; coord(3,4) = 0.0d0
    coord(1,5) = 0.0d0; coord(2,5) = 0.0d0; coord(3,5) = 1.0d0
    coord(1,6) = 1.0d0; coord(2,6) = 0.0d0; coord(3,6) = 1.0d0
    coord(1,7) = 1.0d0; coord(2,7) = 1.0d0; coord(3,7) = 1.0d0
    coord(1,8) = 0.0d0; coord(2,8) = 1.0d0; coord(3,8) = 1.0d0

    !> case 1
    global_pos(1) = 0.0d0
    global_pos(2) = 0.0d0
    global_pos(3) = 0.0d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) =-1.0d0
    ans(2) =-1.0d0
    ans(3) =-1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 1", local_pos, ans)

    !> case 2
    global_pos(1) = 1.0d0
    global_pos(2) = 1.0d0
    global_pos(3) = 1.0d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 1.0d0
    ans(2) = 1.0d0
    ans(3) = 1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 2", local_pos, ans)

    !> case 3
    global_pos(1) = 0.5d0
    global_pos(2) = 0.5d0
    global_pos(3) = 0.5d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 0.0d0
    ans(2) = 0.0d0
    ans(3) = 0.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 3", local_pos, ans)

    !> pattern 2
    coord(1,1) = 0.0d0; coord(2,1) = 0.0d0; coord(3,1) = 0.0d0
    coord(1,2) = 1.0d0; coord(2,2) = 0.0d0; coord(3,2) = 0.0d0
    coord(1,3) = 1.0d0; coord(2,3) = 1.0d0; coord(3,3) = 0.0d0
    coord(1,4) = 0.0d0; coord(2,4) = 1.0d0; coord(3,4) = 0.0d0
    coord(1,5) = 1.0d0; coord(2,5) = 0.0d0; coord(3,5) = 1.0d0
    coord(1,6) = 2.0d0; coord(2,6) = 0.0d0; coord(3,6) = 1.0d0
    coord(1,7) = 1.0d0; coord(2,7) = 1.0d0; coord(3,7) = 1.0d0
    coord(1,8) = 0.0d0; coord(2,8) = 1.0d0; coord(3,8) = 1.0d0

    !> case 4
    global_pos(1) = 0.0d0
    global_pos(2) = 0.0d0
    global_pos(3) = 0.0d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) =-1.0d0
    ans(2) =-1.0d0
    ans(3) =-1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 4", local_pos, ans)

    !> case 5
    global_pos(1) = 1.0d0
    global_pos(2) = 1.0d0
    global_pos(3) = 1.0d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) = 1.0d0
    ans(2) = 1.0d0
    ans(3) = 1.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 5", local_pos, ans)

    !> case 6
    global_pos(1) = 0.5d0
    global_pos(2) = 0.5d0
    global_pos(3) = 0.5d0

    call gg_tools_get_distance_3d(n_base, coord, fptr, &
      & global_pos, local_pos, ths, ths_up, eps, is_converge)

    ans(1) =-0.5d0
    ans(2) = 0.0d0
    ans(3) = 0.0d0

    call monolis_test_check_eq_R("gg_tools_get_distance_3d_test 6", local_pos, ans)
  end subroutine gg_tools_get_distance_3d_test

  subroutine get_shapefunc_deriv_3d_test()
    implicit none
    integer(kint) :: n_base
    real(kdouble) :: local_pos(3), eps
    real(kdouble) :: deriv(8,3), ans(8,3)

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !> 要素を構成する形状関数の個数
        integer(4), intent(in) :: n_base
        !> 入力座標に最も近い局所座標
        real(8), intent(in) :: local_pos(3)
        !> 形状関数の重み
        real(8), intent(out) :: weight(n_base)
      end subroutine shape_func
    end interface

    procedure(shape_func), pointer ::  fptr => null()

    call monolis_std_global_log_string("get_shapefunc_deriv_3d")

    n_base = 8
    eps = 1.0d-3

    fptr => C3D8_shapefunc

    !> case 1
    local_pos(1) = 0.0d0
    local_pos(2) = 0.0d0
    local_pos(3) = 0.0d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 1a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 1b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 1c", deriv(:,3), ans(:,3))

    !> case 2
    local_pos(1) =-0.5d0
    local_pos(2) =-0.5d0
    local_pos(3) =-0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 2a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 2b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 2c", deriv(:,3), ans(:,3))

    !> case 3
    local_pos(1) = 0.5d0
    local_pos(2) =-0.5d0
    local_pos(3) =-0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 3a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 3b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 3c", deriv(:,3), ans(:,3))

    !> case 4
    local_pos(1) = 0.5d0
    local_pos(2) = 0.5d0
    local_pos(3) =-0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 4a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 4b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 4c", deriv(:,3), ans(:,3))

    !> case 5
    local_pos(1) =-0.5d0
    local_pos(2) = 0.5d0
    local_pos(3) =-0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 5a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 5b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 5c", deriv(:,3), ans(:,3))

    !> case 6
    local_pos(1) =-0.5d0
    local_pos(2) =-0.5d0
    local_pos(3) = 0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 6a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 6b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 6c", deriv(:,3), ans(:,3))

    !> case 7
    local_pos(1) = 0.5d0
    local_pos(2) =-0.5d0
    local_pos(3) = 0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 7a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 7b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 7c", deriv(:,3), ans(:,3))

    !> case 8
    local_pos(1) = 0.5d0
    local_pos(2) = 0.5d0
    local_pos(3) = 0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 8a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 8b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 8c", deriv(:,3), ans(:,3))

    !> case 9
    local_pos(1) =-0.5d0
    local_pos(2) = 0.5d0
    local_pos(3) = 0.5d0

    call get_shapefunc_deriv_3d(n_base, fptr, local_pos, eps, deriv)
    call C3D8_shapefunc_deriv(n_base, local_pos, ans)

    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 9a", deriv(:,1), ans(:,1))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 9b", deriv(:,2), ans(:,2))
    call monolis_test_check_eq_R("get_shapefunc_deriv_3d 9c", deriv(:,3), ans(:,3))
  end subroutine get_shapefunc_deriv_3d_test

  subroutine C2D4_shapefunc(n_base, local, func)
    implicit none
    integer(kint), intent(in) :: n_base
    real(kdouble), intent(in) :: local(2)
    real(kdouble), intent(out) :: func(n_base)

    func(1) = 0.25d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(2) = 0.25d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(3) = 0.25d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(4) = 0.25d0*(1.0d0-local(1))*(1.0d0+local(2))
  end subroutine C2D4_shapefunc

  subroutine C2D4_shapefunc_deriv(n_base, local, func)
    implicit none
    integer(kint), intent(in) :: n_base
    real(kdouble), intent(in) :: local(2)
    real(kdouble), intent(out) :: func(n_base,2)

    func(1,1) = -0.25d0*(1.0d0-local(2))
    func(2,1) =  0.25d0*(1.0d0-local(2))
    func(3,1) =  0.25d0*(1.0d0+local(2))
    func(4,1) = -0.25d0*(1.0d0+local(2))

    func(1,2) = -0.25d0*(1.0d0-local(1))
    func(2,2) = -0.25d0*(1.0d0+local(1))
    func(3,2) =  0.25d0*(1.0d0+local(1))
    func(4,2) =  0.25d0*(1.0d0-local(1))
  end subroutine C2D4_shapefunc_deriv

  subroutine C3D8_shapefunc(n_base, local, func)
    implicit none
    integer(kint), intent(in) :: n_base
    real(kdouble), intent(in) :: local(3)
    real(kdouble), intent(out) :: func(n_base)

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
    real(kdouble), intent(out) :: func(n_base,3)

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
end module mod_gg_tools_distance_determination_test
