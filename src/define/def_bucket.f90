!> バケットセルモジュール
module mod_def_bucket
  use mod_monolis_utils
  implicit none

  !> バケットセル構造体
  type type_monolis_neighbor_search_main
    !> バケットセルに登録された番号数
    integer(kint) :: nid = 0
    !> バケットセルに登録された番号
    integer(kint), allocatable :: id(:)
    !> memo: 次元ごとに番号を記憶する
  end type type_monolis_neighbor_search_main

  !> バケット構造体
  type type_monolis_neighbor_search
    !> バケットの分割数（nx, ny, nz）
    integer(kint) :: div(3)
    !> バウンダリボックス（xmin, xmax, ymin, ymax, zmin, zmax）
    real(kdouble) :: BB(6)
    !> セルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    !> バケットセル構造体（配列サイズ nx × ny × nz）
    type(type_monolis_neighbor_search_main), allocatable :: cell(:)
  end type type_monolis_neighbor_search

end module mod_def_bucket
