!> 距離判定モジュール
module mod_gg_tools_distance_determination
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup distance
  !> 距離判定関数
  subroutine gg_tools_get_distance_2d(n_base, coord, shape_func, &
    & global_pos, local_pos, ths, ths_up, eps, is_converge)
    implicit none
    !> 要素を構成する形状関数の個数
    integer(kint), intent(in) :: n_base
    !> 要素を構成する座標
    real(kdouble), intent(in) :: coord(3,n_base)
    !> 形状関数の関数ポインタ
    pointer :: shape_func
    !> 距離判定したい入力座標
    real(kdouble), intent(in) :: global_pos(3)
    !> 入力座標に最も近い局所座標
    real(kdouble), intent(out) :: local_pos(2)
    !> 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> 発散判定閾値
    real(kdouble), intent(in) :: ths_up
    !> 数値微分の差分値
    real(kdouble), intent(in) :: eps
    !> 収束判定フラグ
    logical, intent(out) :: is_converge

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

    integer(kint) :: i, j
    real(kdouble) :: tangent(3,2), invJacob(2,2), det
    real(kdouble) :: norm, n(n_base), deriv(n_base,2), r(3), dr(2), dF(2), d2F(2,2)
    logical :: is_fail

    is_converge = .false.
    local_pos = 0.0d0
    r = 0.0d0

    do i = 1, 10
      call shape_func(n_base, local_pos, n)
      r = matmul(coord, n) - global_pos
      norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
      if(norm < ths)then
        is_converge = .true.
        exit
      endif

      call get_shapefunc_deriv_2d(n_base, shape_func, local_pos, eps, deriv)

      tangent = matmul(coord, deriv)
      dF = - matmul(r, tangent)
      d2F(1,1) = dot_product(tangent(:,1), tangent(:,1))
      d2F(1,2) = dot_product(tangent(:,1), tangent(:,2))
      d2F(2,1) = dot_product(tangent(:,2), tangent(:,1))
      d2F(2,2) = dot_product(tangent(:,2), tangent(:,2))

      call monolis_get_inverse_matrix_R_2d(d2F, invJacob, det, is_fail)

      if(is_fail) exit
      dr = matmul(invJacob, dF)
      local_pos = local_pos + dr

      norm = dsqrt(local_pos(1)**2 + local_pos(2)**2)
      if(ths_up < norm)then
        exit
      endif
    enddo
  end subroutine gg_tools_get_distance_2d

  !> @ingroup distance
  !> 距離判定関数
  subroutine gg_tools_get_distance_3d(n_base, coord, shape_func, &
    & global_pos, local_pos, ths, ths_up, eps, is_converge)
    implicit none
    !> 要素を構成する形状関数の個数
    integer(kint), intent(in) :: n_base
    !> 要素を構成する座標
    real(kdouble), intent(in) :: coord(3,n_base)
    !> 形状関数の関数ポインタ
    pointer :: shape_func
    !> 距離判定したい入力座標
    real(kdouble), intent(in) :: global_pos(3)
    !> 入力座標に最も近い局所座標
    real(kdouble), intent(out) :: local_pos(3)
    !> 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> 発散判定閾値
    real(kdouble), intent(in) :: ths_up
    !> 数値微分の差分値
    real(kdouble), intent(in) :: eps
    !> 収束判定フラグ
    logical, intent(out) :: is_converge

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

    integer(kint) :: i, j
    real(kdouble) :: jacobi(3,3), invJacob(3,3), det
    real(kdouble) :: norm, n(n_base), deriv(n_base,3), r(3), dr(3)
    logical :: is_fail

    is_converge = .false.
    local_pos = 0.0d0
    r = 0.0d0

    do i = 1, 10
      call shape_func(n_base, local_pos, n)
      r = matmul(coord, n) - global_pos
      norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
      if(norm < ths)then
        is_converge = .true.
        exit
      endif

      call get_shapefunc_deriv_3d(n_base, shape_func, local_pos, eps, deriv)

      jacobi = matmul(coord, deriv)
      call monolis_get_inverse_matrix_R_3d(jacobi, invJacob, det, is_fail)
      if(is_fail) exit

      dr = - matmul(invJacob, r)
      local_pos = local_pos + dr

      norm = dsqrt(local_pos(1)**2 + local_pos(2)**2 + local_pos(3)**2)
      if(ths_up < norm)then
        exit
      endif
    enddo
  end subroutine gg_tools_get_distance_3d

  !> @ingroup distance
  !> 形状関数の微分値の取得（中心差分）
  subroutine get_shapefunc_deriv_3d(n_base, shape_func, local_pos, eps, deriv)
    implicit none
    !> 要素を構成する形状関数の個数
    integer(kint), intent(in) :: n_base
    !> 形状関数の関数ポインタ
    pointer :: shape_func
    !> 入力座標に最も近い局所座標
    real(kdouble), intent(in) :: local_pos(3)
    !> 数値微分の差分値
    real(kdouble), intent(in) :: eps
    !> 形状関数の微分
    real(kdouble), intent(out) :: deriv(n_base,3)

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

    integer(kint) :: dim
    real(kdouble) :: pos(3)
    real(kdouble) :: n1(n_base), n2(n_base)

    do dim = 1, 3
      pos = local_pos
      pos(dim) = pos(dim) - eps
      call shape_func(n_base, pos, n1)

      pos = local_pos
      pos(dim) = pos(dim) + eps
      call shape_func(n_base, pos, n2)

      deriv(:,dim) = (n2 - n1)/(2.0d0*eps)
    enddo
  end subroutine get_shapefunc_deriv_3d

  !> @ingroup distance
  !> 形状関数の微分値の取得（中心差分）
  subroutine get_shapefunc_deriv_2d(n_base, shape_func, local_pos, eps, deriv)
    implicit none
    !> 要素を構成する形状関数の個数
    integer(kint), intent(in) :: n_base
    !> 形状関数の関数ポインタ
    pointer :: shape_func
    !> 入力座標に最も近い局所座標
    real(kdouble), intent(in) :: local_pos(2)
    !> 数値微分の差分値
    real(kdouble), intent(in) :: eps
    !> 形状関数の微分
    real(kdouble), intent(out) :: deriv(n_base,2)

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

    integer(kint) :: dim
    real(kdouble) :: pos(2)
    real(kdouble) :: n1(n_base), n2(n_base)

    do dim = 1, 2
      pos = local_pos
      pos(dim) = pos(dim) - eps
      call shape_func(n_base, pos, n1)

      pos = local_pos
      pos(dim) = pos(dim) + eps
      call shape_func(n_base, pos, n2)

      deriv(:,dim) = (n2 - n1)/(2.0d0*eps)
    enddo
  end subroutine get_shapefunc_deriv_2d
end module mod_gg_tools_distance_determination
