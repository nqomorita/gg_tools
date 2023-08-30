!> バケットセルモジュール
module mod_gg_tools_def_bucket
  use mod_monolis_utils
  implicit none

  !> @ingroup bucket
  !> バケットセル構造体
  type type_gg_tools_bucket_search_main
    !> バケットセルに登録された番号数
    integer(kint) :: nid = 0
    !> バケットセルに登録された番号
    integer(kint), allocatable :: id(:)
  end type type_gg_tools_bucket_search_main

  !> @ingroup bucket
  !> バケット検索構造体
  type type_gg_tools_bucket_search
    !> バケットの分割数（nx, ny, nz）
    integer(kint) :: n_div(3)
    !> バウンダリボックス（xmin, xmax, ymin, ymax, zmin, zmax）
    real(kdouble) :: BB(6)
    !> セルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    !> 検索領域の判定閾値
    real(kdouble) :: ths
    !> 要素次元ごとのバケットセル構造体（配列サイズ nx × ny × nz）
    type(type_gg_tools_bucket_search_main), allocatable :: cell_0d(:)
    type(type_gg_tools_bucket_search_main), allocatable :: cell_1d(:)
    type(type_gg_tools_bucket_search_main), allocatable :: cell_2d(:)
    type(type_gg_tools_bucket_search_main), allocatable :: cell_3d(:)
  end type type_gg_tools_bucket_search

contains

  !> @ingroup bucket
  !> バケット検索の初期化
  subroutine gg_tools_bucket_search_init(ggt_bucket_search, BB, n_div)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 検索領域を定義するバウンダリボックス（xmin, xmax, ymin, ymax, zmin, zmax）
    real(kdouble) :: BB(6)
    !> バケットの分割数（nx, ny, nz）
    integer(kint) :: n_div(3)
    integer(kint) :: n_all
    real(kdouble) :: minbb

    ggt_bucket_search%n_div = n_div
    ggt_bucket_search%BB = BB
    ggt_bucket_search%dx(1) = (BB(2)-BB(1))/dble(n_div(1))
    ggt_bucket_search%dx(2) = (BB(4)-BB(3))/dble(n_div(2))
    ggt_bucket_search%dx(3) = (BB(6)-BB(5))/dble(n_div(3))

    minbb = min(BB(2)-BB(1), BB(4)-BB(3), BB(6)-BB(5))
    ggt_bucket_search%ths = 1.0d-3*minbb

    n_all = n_div(1)*n_div(2)*n_div(3)
    allocate(ggt_bucket_search%cell_0d(n_all))
    allocate(ggt_bucket_search%cell_1d(n_all))
    allocate(ggt_bucket_search%cell_2d(n_all))
    allocate(ggt_bucket_search%cell_3d(n_all))
  end subroutine gg_tools_bucket_search_init

  !> @ingroup bucket
  !> バケットへの情報登録
  subroutine gg_tools_bucket_search_push(ggt_bucket_search, BB, id)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 要素領域を定義するバウンダリボックス（xmin, xmax, ymin, ymax, zmin, zmax）
    real(kdouble) :: BB(6)
    !> 要素領域 id
    integer(kint) :: id
    integer(kint) :: x, y, z, in, imin(3), imax(3)
    real(kdouble) :: pos(3), ths

    ths = ggt_bucket_search%ths

    pos(1) = BB(1) - ths
    pos(2) = BB(3) - ths
    pos(3) = BB(5) - ths
    call get_int_coordinate(ggt_bucket_search, pos, imin)

    pos(1) = BB(2) + ths
    pos(2) = BB(4) + ths
    pos(3) = BB(6) + ths
    call get_int_coordinate(ggt_bucket_search, pos, imax)

    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      in = get_index(ggt_bucket_search%n_div, x, y, z)
      call gg_tools_bucket_search_push_main(ggt_bucket_search, in, id)
    enddo
    enddo
    enddo
  end subroutine gg_tools_bucket_search_push

  !> @ingroup dev
  !> バケットへの情報登録（メイン関数）
  subroutine gg_tools_bucket_search_push_main(ggt_bucket_search, in, id)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> バケットセルの id
    integer(kint) :: in
    !> 登録する要素領域 id
    integer(kint) :: id
    integer(kint) :: add(1)

    add = id
    call monolis_append_I_1d(ggt_bucket_search%cell_3d(in)%id, 1, add)
    ggt_bucket_search%cell_3d(in)%nid = ggt_bucket_search%cell_3d(in)%nid + 1
  end subroutine gg_tools_bucket_search_push_main

  !> @ingroup bucket
  !> バケットから登録情報の取得（座標を入力）
  subroutine gg_tools_bucket_search_get_by_position(ggt_bucket_search, pos, nid, id)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 検索座標
    real(kdouble) :: pos(3)
    !> 取得情報の個数
    integer(kint) :: nid
    !> 取得情報を保持した整数配列
    integer(kint), allocatable :: id(:)
    integer(kint) :: in
    logical :: is_in

    nid = 0
    call monolis_dealloc_I_1d(id)

    call BB_check(ggt_bucket_search, pos, is_in)
    if(.not. is_in) return

    call get_int_coordinate(ggt_bucket_search, pos, id)
    in = get_index(ggt_bucket_search%n_div, id(1), id(2), id(3))

    nid = ggt_bucket_search%cell_3d(in)%nid

    if(nid > 0)then
      call monolis_alloc_I_1d(id, nid)
      id = ggt_bucket_search%cell_3d(in)%id
    endif
  end subroutine gg_tools_bucket_search_get_by_position

  !> @ingroup bucket
  !> バケットから登録情報の取得（バウンダリボックスを入力）
  subroutine gg_tools_bucket_search_get_by_bb(ggt_bucket_search, BB, nid, id)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 検索領域を定義するバウンダリボックス（xmin, xmax, ymin, ymax, zmin, zmax）
    real(kdouble) :: BB(6)
    !> 取得情報の個数
    integer(kint) :: nid
    !> 取得情報を保持した整数配列
    integer(kint), allocatable :: id(:)
    integer(kint) :: imin(3), imax(3)
    integer(kint) :: in, x, y, z, newlen
    integer(kint), allocatable :: tmp(:)
    real(kdouble) :: pos(3), ths
    logical :: is_in

    nid = 0
    call monolis_dealloc_I_1d(id)

    call BB_modify(ggt_bucket_search, BB, is_in)
    if(.not. is_in) return

    ths = ggt_bucket_search%ths

    pos(1) = BB(1) - ths
    pos(2) = BB(3) - ths
    pos(3) = BB(5) - ths
    call get_int_coordinate(ggt_bucket_search, pos, imin)

    pos(1) = BB(2) + ths
    pos(2) = BB(4) + ths
    pos(3) = BB(6) + ths
    call get_int_coordinate(ggt_bucket_search, pos, imax)

    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      in = get_index(ggt_bucket_search%n_div, x, y, z)
      if(ggt_bucket_search%cell_3d(in)%nid > 0)then
        call monolis_append_I_1d(tmp, ggt_bucket_search%cell_3d(in)%nid, ggt_bucket_search%cell_3d(in)%id)
        nid = nid + ggt_bucket_search%cell_3d(in)%nid
      endif
    enddo
    enddo
    enddo

    if(.not. allocated(tmp)) return

    call monolis_qsort_I_1d(tmp, 1, nid)
    call monolis_get_uniq_array_I(tmp, nid, newlen)

    call monolis_alloc_I_1d(id, newlen)
    id = tmp(1:newlen)
    call monolis_dealloc_I_1d(tmp)
    nid = newlen
  end subroutine gg_tools_bucket_search_get_by_bb

  !> @ingroup bucket
  !> バケット情報の終了処理
  subroutine gg_tools_bucket_search_finalize(ggt_bucket_search)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    deallocate(ggt_bucket_search%cell_0d)
    deallocate(ggt_bucket_search%cell_1d)
    deallocate(ggt_bucket_search%cell_2d)
    deallocate(ggt_bucket_search%cell_3d)
  end subroutine gg_tools_bucket_search_finalize

  !> @ingroup dev
  !> 入力座標とバウンディングボックスの内包判定
  subroutine BB_check(ggt_bucket_search, pos, is_in)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 検索座標
    real(kdouble) :: pos(3)
    !> 内包フラグ
    logical :: is_in
    real(kdouble) :: BB(6), ths

    is_in = .false.
    BB = ggt_bucket_search%BB
    ths = ggt_bucket_search%ths

    if( BB(1)-ths <= pos(1) .and. pos(1) <= BB(2)+ths .and. &
      & BB(3)-ths <= pos(2) .and. pos(2) <= BB(4)+ths .and. &
      & BB(5)-ths <= pos(3) .and. pos(3) <= BB(6)+ths )then
      is_in = .true.
    endif
  end subroutine BB_check

  !> @ingroup dev
  !> 検索用バウンディングボックスの範囲修正
  subroutine BB_modify(ggt_bucket_search, BB_in, is_in)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 入力バウンディングボックス
    real(kdouble) :: BB_in(6)
    !> 内包フラグ
    logical :: is_in
    real(kdouble) :: BB(6)

    is_in = .true.
    BB = ggt_bucket_search%BB

    if(BB_in(2) < BB(1)) is_in = .false.
    if(BB_in(4) < BB(3)) is_in = .false.
    if(BB_in(6) < BB(5)) is_in = .false.
    if(BB_in(1) > BB(2)) is_in = .false.
    if(BB_in(3) > BB(4)) is_in = .false.
    if(BB_in(5) > BB(6)) is_in = .false.

    if(.not. is_in) return

    if(BB_in(1) < BB(1)) BB_in(1) = BB(1)
    if(BB_in(3) < BB(3)) BB_in(3) = BB(3)
    if(BB_in(5) < BB(5)) BB_in(5) = BB(5)
    if(BB_in(2) > BB(2)) BB_in(2) = BB(2)
    if(BB_in(4) > BB(4)) BB_in(4) = BB(4)
    if(BB_in(6) > BB(6)) BB_in(6) = BB(6)
  end subroutine BB_modify

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

  !> @ingroup dev
  !> バケットセル整数座標の取得
  subroutine get_int_coordinate(ggt_bucket_search, pos, id)
    implicit none
    !> バケット検索構造体
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    !> 入力座標
    real(kdouble) :: pos(3)
    !> バケットセル整数座標
    integer(kint) :: id(3)

    id(1) = int((pos(1) - ggt_bucket_search%BB(1))/ggt_bucket_search%dx(1)) + 1
    id(2) = int((pos(2) - ggt_bucket_search%BB(3))/ggt_bucket_search%dx(2)) + 1
    id(3) = int((pos(3) - ggt_bucket_search%BB(5))/ggt_bucket_search%dx(3)) + 1
    if(id(1) < 1) id(1) = 1
    if(id(2) < 1) id(2) = 1
    if(id(3) < 1) id(3) = 1
    if(id(1) > ggt_bucket_search%n_div(1)) id(1) = ggt_bucket_search%n_div(1)
    if(id(2) > ggt_bucket_search%n_div(2)) id(2) = ggt_bucket_search%n_div(2)
    if(id(3) > ggt_bucket_search%n_div(3)) id(3) = ggt_bucket_search%n_div(3)
  end subroutine get_int_coordinate
end module mod_gg_tools_def_bucket
