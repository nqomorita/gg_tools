!> バケットセルモジュール
module mod_ggtools_def_bucket
  use mod_monolis_utils
  implicit none

  real(kdouble) :: ths = 1.0d-6

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
  !> @details バウンダリボックスの最小座標 xmin を起点に、バウンディングボックスを覆う最小のバケットセルを確保する。
  !> @details バケットセル確保後に、バウンダリボックスの最大座標 xmax が更新される。
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

    ggtools_bucket%xmin = xmin
    ggtools_bucket%xmax = xmax
    ggtools_bucket%dx = dx

    ggtools_bucket%nx(1) = 1
    ggtools_bucket%nx(2) = 2
    ggtools_bucket%nx(3) = 3
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

    ggtools_bucket%xmin = xmin
    ggtools_bucket%xmax = xmax
    ggtools_bucket%nx = nx

    ggtools_bucket%dx(1) = 1.0d0
    ggtools_bucket%dx(2) = 2.0d0
    ggtools_bucket%dx(3) = 3.0d0
  end subroutine ggtools_bucket_init_by_nx

  !> @ingroup bucket
  !> バケットの終了処理
  subroutine ggtools_bucket_finalize(ggtools_bucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket

    ggtools_bucket%xmin = 0.0d0
    ggtools_bucket%xmax = 0.0d0
    ggtools_bucket%dx = 0.0d0
    ggtools_bucket%nx = 0
  end subroutine ggtools_bucket_finalize

  !> @ingroup bucket
  !> バケットセルの初期化処理
  subroutine ggtools_bucket_cell_init(ggtools_bucket_cell, nx)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_bucket_cell), allocatable :: ggtools_bucket_cell(:)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    integer(kint) :: n_total

    n_total = nx(1)*nx(2)*nx(3)
    allocate(ggtools_bucket_cell(n_total))
  end subroutine ggtools_bucket_cell_init

  !> @ingroup bucket
  !> バケットセルの終了処理
  subroutine ggtools_bucket_cell_finalize(ggtools_bucket_cell)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell(:)
    integer(kint) :: i

    do i = 1, size(ggtools_bucket_cell)
      if(ggtools_bucket_cell(i)%nid == 0) cycle
      deallocate(ggtools_bucket_cell(i)%id)
    enddo
  end subroutine ggtools_bucket_cell_finalize

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
    integer(kint) :: x, y, z, in, imin(3), imax(3)
    real(kdouble) :: x_point(3)

    x_point(1) = xmin_local(1) - ths
    x_point(2) = xmin_local(2) - ths
    x_point(3) = xmin_local(3) - ths
    call ggtools_bucket_get_int_coordinate(ggtools_bucket, x_point, imin)

    x_point(1) = xmax_local(1) + ths
    x_point(2) = xmax_local(2) + ths
    x_point(3) = xmax_local(3) + ths
    call ggtools_bucket_get_int_coordinate(ggtools_bucket, x_point, imax)

    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      in = get_index(ggtools_bucket%nx, x, y, z)
      call ggtools_bucket_set_id_main(ggtools_bucket_cell(in), id)
    enddo
    enddo
    enddo
  end subroutine ggtools_bucket_set_id_by_lbb

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
    integer(kint) :: int_id(3)

    call get_int_coordinate(ggtools_bucket, x_point, int_id)
  end subroutine ggtools_bucket_set_id_by_point

  !> @ingroup dev
  !> バケットへの情報登録（メイン関数）
  subroutine ggtools_bucket_set_id_main(ggtools_bucket_cell, data)
    implicit none
    !> バケット検索構造体
    type(type_ggtools_bucket_cell) :: ggtools_bucket_cell
    !> 登録する要素領域 id
    integer(kint) :: data
    integer(kint) :: add(1)

    add = data
    call monolis_append_I_1d(ggtools_bucket_cell%id, 1, add)
    ggtools_bucket_cell%nid = ggtools_bucket_cell%nid + 1
  end subroutine ggtools_bucket_set_id_main

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
  end subroutine ggtools_bucket_get_by_point

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
  end subroutine ggtools_bucket_get_by_lbb

  !> @ingroup dev
  !> バケットセル整数座標の取得
  subroutine ggtools_bucket_get_int_coordinate(ggtools_bucket, pos, id)
    implicit none
    !> バケット検索構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> 入力座標
    real(kdouble) :: pos(3)
    !> バケットセル整数座標
    integer(kint) :: id(3)

    id(1) = int((pos(1) - ggtools_bucket%xmin(1))/ggtools_bucket%dx(1)) + 1
    id(2) = int((pos(2) - ggtools_bucket%xmin(2))/ggtools_bucket%dx(2)) + 1
    id(3) = int((pos(3) - ggtools_bucket%xmin(3))/ggtools_bucket%dx(3)) + 1
    if(id(1) < 1) id(1) = 1
    if(id(2) < 1) id(2) = 1
    if(id(3) < 1) id(3) = 1
    if(id(1) > ggtools_bucket%nx(1)) id(1) = ggtools_bucket%nx(1)
    if(id(2) > ggtools_bucket%nx(2)) id(2) = ggtools_bucket%nx(2)
    if(id(3) > ggtools_bucket%nx(3)) id(3) = ggtools_bucket%nx(3)
  end subroutine ggtools_bucket_get_int_coordinate

  !> @ingroup dev
  !> バケットセル整数座標からバケットセル id を取得
  function get_index(div, x, y, z)
    implicit none
    !> バケットセル分割数
    integer(kint) :: div(3)
    !> バケットセル整数座標 x
    integer(kint) :: x
    !> バケットセル整数座標 y
    integer(kint) :: y
    !> バケットセル整数座標 z
    integer(kint) :: z
    !> バケットセル id
    integer(kint) :: get_index
    get_index = x + (y-1)*div(1) + (z-1)*div(1)*div(2)
  end function get_index

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
    xmin = ggtools_bucket%xmin
    xmax = ggtools_bucket%xmax
  end subroutine ggtools_bucket_get_bucket_size

  !> @ingroup bucket
  !> バケットセル分割数を取得
  subroutine ggtools_bucket_get_number_of_bucket_divisions(ggtools_bucket, nx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    nx = ggtools_bucket%nx
  end subroutine ggtools_bucket_get_number_of_bucket_divisions

  !> @ingroup bucket
  !> バケットセルサイズを取得
  subroutine ggtools_bucket_get_bucket_cell_size(ggtools_bucket, dx)
    implicit none
    !> バケット構造体
    type(type_ggtools_bucket) :: ggtools_bucket
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    dx = ggtools_bucket%dx
  end subroutine ggtools_bucket_get_bucket_cell_size
end module mod_ggtools_def_bucket
