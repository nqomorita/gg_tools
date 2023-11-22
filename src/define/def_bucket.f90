!> バケットセルモジュール
module mod_ggtools_def_bucket
  use mod_monolis_utils
  implicit none

  !> @ingroup bucket
  !> 第 i 番目のバケットセルの基本情報構造体
  type type_ggtools_bucket_cell
    !> バケットセルに登録された個数
    integer(kint) :: nid = 0
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id(:)
  end type type_ggtools_bucket_cell

  !> @ingroup bucket
  !> バケットの基本情報構造体
  type type_ggtools_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
  end type type_ggtools_bucket

contains

  !> @ingroup bucket
  !> バケットの初期化処理（バケットセルサイズを入力）
  subroutine ggtools_bucket_init(ggtools_bucket, xmin, xmax, dx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
  end subroutine ggtools_bucket_init

  !> @ingroup bucket
  !> バケットの初期化処理（バケットセル分割数を入力）
  subroutine ggtools_bucket_init_by_nx(ggtools_bucket, xmin, xmax, nx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
  end subroutine ggtools_bucket_init_by_nx

  !> @ingroup bucket
  !> バケットの終了処理
  subroutine ggtools_bucket_finalize(ggtools_bucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
  end subroutine

  !> @ingroup bucket
  !> バケットセルの初期化処理
  subroutine ggtools_bucket_cell_init(ggtools_bucket_cell, nx)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
  end subroutine

  !> @ingroup bucket
  !> バケットセルの終了処理
  subroutine ggtools_bucket_cell_finalize(ggtools_bucket_cell)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
  end subroutine

  !> @ingroup bucket
  !> LBB を指定してバケットセルに id を登録する関数
  subroutine ggtools_bucket_set_id_by_lbb(ggtools_bucket, ggtools_bucket_cell, xmin_local, xmax_local, id)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    !> ローカルバウンダリボックス LBB の最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンダリボックス LBB の最大座標
    real(kdouble) :: xmax_local(3)
    !> バケットセルに登録する id
    integer(kint) :: id
  end subroutine

  !> @ingroup bucket
  !> 座標を指定して 1 つのバケットセルに id を登録する関数
  subroutine ggtools_bucket_set_id_by_point(ggtools_bucket, ggtools_bucket_cell, x_point, id)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バケットセルに登録する id
    integer(kint) :: id
  end subroutine

  !> @ingroup bucket
  !> 入力した座標を内包するバケットセルから、nid 個の整数配列 id_array を取得する関数
  subroutine ggtools_bucket_get_by_point(ggtools_bucket, ggtools_bucket_cell, x_point, nid, id_array)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
  end subroutine

  !> @ingroup bucket
  !> 入力した LBB を内包するバケットセルから、nid 個の整数配列 id_array を取得する関数
  subroutine ggtools_bucket_get_by_lbb(ggtools_bucket, ggtools_bucket_cell, xmin_local, xmax_local, nid, id_array)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    !> ローカルバウンダリボックス LBB の最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンダリボックス LBB の最大座標
    real(kdouble) :: xmax_local(3)
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
    !内側の処理
    !lbb から lbb を含む整数座標の int の min, max を取得（x, y, z の 3 方向）
    !バケットセルの番号列を返す
    !int の min, max の範囲で、バケットセルにアクセスし、個々の n_id_i, id_array_i を取得（ただし実体としては持っていない）
    !n_id_i, id_array_i から重複削除し、n_id, id_array を計算
  end subroutine

  !> @ingroup bucket
  !> バケットサイズを取得
  subroutine ggtools_bucket_get_bucket_size(ggtools_bucket, xmin, xmax)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
  end subroutine ggtools_bucket_get_bucket_size

  !> @ingroup bucket
  !> バケットセル分割数を取得
  subroutine ggtools_bucket_get_number_of_bucket_divisions(ggtools_bucket, nx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
  end subroutine ggtools_bucket_get_number_of_bucket_divisions

  !> @ingroup bucket
  !> バケットセルサイズを取得
  subroutine ggtools_bucket_get_bucket_cell_size(ggtools_bucket, dx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
  end subroutine ggtools_bucket_get_bucket_cell_size
end module mod_ggtools_def_bucket
