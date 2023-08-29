!> 距離判定モジュール
module mod_gg_tools_distance_determination
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup distance
  !> 距離判定関数
  subroutine gg_tools_get_distance_3d(n_base, coord, shape_func, &
    & global_pos, local_pos, ths, ths_up, is_converge)
    implicit none
    !>
    integer(kint), intent(in) :: n_base
    !>
    real(kdouble), intent(in) :: coord(3,n_base)
    !>
    pointer :: shape_func
    !>
    real(kdouble), intent(in) :: global_pos(3)
    !>
    real(kdouble), intent(out) :: local_pos(3)
    !>
    real(kdouble), intent(in) :: ths
    !>
    real(kdouble), intent(in) :: ths_up
    !>
    logical, intent(out) :: is_converge
    integer(kint) :: i, j
    real(kdouble) :: jacobi(3,3), invJacob(3,3), det
    real(kdouble) :: norm, n(n_base), deriv(n_base,3), r(3), dr(3)
    logical :: is_fail

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !>
        integer(4), intent(in) :: n_base
        !>
        real(8), intent(in) :: local_pos(3)
        !>
        real(8), intent(out) :: weight(n_base, 3)
      end subroutine shape_func
    end interface

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

      call get_shapefunc_deriv(n_base, shape_func, local_pos, deriv)

      jacobi = matmul(coord, deriv)
      call monolis_get_inverse_matrix_R_3d(jacobi, invJacob, det, is_fail)
      if(is_fail) exit

      dr = - Matmul(invJacob, r)
      do j = 1, 3
        local_pos(j) = local_pos(j) + dr(j)
      enddo

      norm = dsqrt(local_pos(1)**2 + local_pos(2)**2 + local_pos(3)**2)
      if(ths_up < norm)then
        exit
      endif
    enddo
  end subroutine gg_tools_get_distance_3d

  !> @ingroup distance
  !> 距離判定関数
  subroutine get_shapefunc_deriv(n_base, shape_func, local_pos, deriv)
    implicit none
    !>
    integer(kint), intent(in) :: n_base
    !>
    pointer :: shape_func
    !>
    real(kdouble), intent(in) :: local_pos(3)
    !>
    real(kdouble), intent(out) :: deriv(n_base,3)

    interface
      subroutine shape_func(n_base, local_pos, weight)
        !>
        integer(4), intent(in) :: n_base
        !>
        real(8), intent(in) :: local_pos(3)
        !>
        real(8), intent(out) :: weight(n_base, 3)
      end subroutine shape_func
    end interface

  end subroutine get_shapefunc_deriv
end module mod_gg_tools_distance_determination
