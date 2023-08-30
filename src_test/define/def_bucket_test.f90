!> バケットセルテストモジュール
module mod_gg_tools_def_bucket_test
  use mod_monolis_utils
  use mod_gg_tools
  implicit none

contains

  subroutine gg_tools_def_bucket_test()
    implicit none

    !# バケット検索用テスト
    call gg_tools_bucket_initialize_test()
    call gg_tools_bucket_search_test()

    !# 一時的なテスト
    call monolis_utils_aabb_test()
    call monolis_utils_kdtree_test()
  end subroutine gg_tools_def_bucket_test

  subroutine gg_tools_bucket_initialize_test()
    implicit none
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    real(kdouble) :: BB(6)
    integer(kint) :: n_div(3)

    call monolis_std_log_string("gg_tools_bucket_search_init")
    call monolis_std_log_string("gg_tools_bucket_search_finalize")

    BB(1) = 0.0d0
    BB(2) = 1.0d0
    BB(3) = 0.0d0
    BB(4) = 1.0d0
    BB(5) = 0.0d0
    BB(6) = 1.0d0

    n_div(1) = 10
    n_div(2) = 10
    n_div(3) = 10

    call gg_tools_bucket_search_init(ggt_bucket_search, BB, n_div)
    call gg_tools_bucket_search_finalize(ggt_bucket_search)
  end subroutine gg_tools_bucket_initialize_test

  subroutine gg_tools_bucket_search_test()
    implicit none
    type(type_gg_tools_bucket_search) :: ggt_bucket_search
    integer(kint) :: n_div(3), id, nid
    real(kdouble) :: BB(6)
    real(kdouble) :: pos(3)
    integer(kint), allocatable :: id_get(:)

    call monolis_std_log_string("gg_tools_bucket_search_push")
    call monolis_std_log_string("gg_tools_bucket_search_get_by_position")

    BB(1) = 0.0d0
    BB(2) = 1.0d0
    BB(3) = 0.0d0
    BB(4) = 1.0d0
    BB(5) = 0.0d0
    BB(6) = 1.0d0

    n_div(1) = 3
    n_div(2) = 3
    n_div(3) = 3

    call gg_tools_bucket_search_init(ggt_bucket_search, BB, n_div)

    BB(1) = 0.0d0
    BB(2) = 0.5d0
    BB(3) = 0.0d0
    BB(4) = 0.5d0
    BB(5) = 0.0d0
    BB(6) = 0.5d0
    id = 1
    call gg_tools_bucket_search_push(ggt_bucket_search, BB, id)

    BB(1) = 0.0d0
    BB(2) = 0.1d0
    BB(3) = 0.0d0
    BB(4) = 0.1d0
    BB(5) = 0.0d0
    BB(6) = 0.1d0
    id = 2
    call gg_tools_bucket_search_push(ggt_bucket_search, BB, id)

    pos(1) = 0.2d0
    pos(2) = 0.2d0
    pos(3) = 0.2d0
    call gg_tools_bucket_search_get_by_position(ggt_bucket_search, pos, nid, id_get)

    call monolis_test_check_eq_I1("gg_tools_bucket_search_test 1a", nid, 2)
    call monolis_test_check_eq_I1("gg_tools_bucket_search_test 1b", id_get(1), 1)
    call monolis_test_check_eq_I1("gg_tools_bucket_search_test 1c", id_get(2), 2)

    pos(1) = 0.5d0
    pos(2) = 0.5d0
    pos(3) = 0.5d0
    call gg_tools_bucket_search_get_by_position(ggt_bucket_search, pos, nid, id_get)

    call monolis_test_check_eq_I1("gg_tools_bucket_search_test 2a", nid, 1)
    call monolis_test_check_eq_I1("gg_tools_bucket_search_test 2b", id_get(1), 1)

    call gg_tools_bucket_search_finalize(ggt_bucket_search)
  end subroutine gg_tools_bucket_search_test

  !> temporary monolis_utils test
  subroutine monolis_utils_aabb_test()
    implicit none

    call monolis_get_aabb_from_coodinates_test()
    call monolis_get_aabb_from_BB_test()
    call monolis_check_inside_in_aabb_test()
  end subroutine monolis_utils_aabb_test

  !> unit test
  subroutine monolis_get_aabb_from_coodinates_test()
    implicit none
    real(kdouble) :: coord(3,3)
    real(kdouble) :: BB(6), BB_ans(6)

    call monolis_std_global_log_string("monolis_get_aabb_from_coodinates")

    coord(1,1) = 1.0d0; coord(2,1) = 10.0d0; coord(3,1) = 100.0d0;
    coord(1,2) = 2.0d0; coord(2,2) = 20.0d0; coord(3,2) = 200.0d0;
    coord(1,3) = 3.0d0; coord(2,3) = 30.0d0; coord(3,3) = 300.0d0;

    call monolis_get_aabb_from_coodinates(coord, BB)

    BB_ans(1) = 1.0d0
    BB_ans(2) = 3.0d0
    BB_ans(3) = 10.0d0
    BB_ans(4) = 30.0d0
    BB_ans(5) = 100.0d0
    BB_ans(6) = 300.0d0

    call monolis_test_check_eq_R("monolis_get_aabb_from_coodinates_test", BB, BB_ans)
  end subroutine monolis_get_aabb_from_coodinates_test

  subroutine monolis_get_aabb_from_BB_test()
    implicit none
    real(kdouble) :: BB_in(6,3)
    real(kdouble) :: BB(6), BB_ans(6)

    call monolis_std_global_log_string("monolis_get_aabb_from_BB")

    BB_in(1,1) = 1.0d0; BB_in(2,1) = 10.0d0; BB_in(3,1) = 4.0d0; BB_in(4,1) = 40.0d0; BB_in(5,1) = 7.0d0; BB_in(6,1) = 70.0d0;
    BB_in(1,2) = 2.0d0; BB_in(2,2) = 20.0d0; BB_in(3,2) = 5.0d0; BB_in(4,2) = 50.0d0; BB_in(5,2) = 8.0d0; BB_in(6,2) = 80.0d0;
    BB_in(1,3) = 3.0d0; BB_in(2,3) = 30.0d0; BB_in(3,3) = 6.0d0; BB_in(4,3) = 60.0d0; BB_in(5,3) = 9.0d0; BB_in(6,3) = 90.0d0;

    call monolis_get_aabb_from_BB(BB_in, BB)

    BB_ans(1) = 1.0d0
    BB_ans(2) = 30.0d0
    BB_ans(3) = 4.0d0
    BB_ans(4) = 60.0d0
    BB_ans(5) = 7.0d0
    BB_ans(6) = 90.0d0

    call monolis_test_check_eq_R("monolis_get_aabb_from_BB_test", BB, BB_ans)
  end subroutine monolis_get_aabb_from_BB_test

  subroutine monolis_check_inside_in_aabb_test()
    implicit none
    real(kdouble) :: coord(3)
    real(kdouble) :: BB(6), ths
    logical :: is_inside

    call monolis_std_global_log_string("monolis_check_inside_in_aabb")

    BB(1) = 1.0d0
    BB(2) = 3.0d0
    BB(3) = 10.0d0
    BB(4) = 30.0d0
    BB(5) = 100.0d0
    BB(6) = 300.0d0

    ths = 0.0d0

    coord(1) = 2.0d0
    coord(2) = 20.0d0
    coord(3) = 200.0d0

    call monolis_check_inside_in_aabb(coord, BB, ths, is_inside)

    call monolis_test_check_eq_L1("monolis_check_inside_in_aabb 1", is_inside, .true.)

    coord(1) = 0.0d0
    coord(2) = 20.0d0
    coord(3) = 200.0d0

    call monolis_check_inside_in_aabb(coord, BB, ths, is_inside)

    call monolis_test_check_eq_L1("monolis_check_inside_in_aabb 2", is_inside, .false.)

    coord(1) = 2.0d0
    coord(2) = 0.0d0
    coord(3) = 200.0d0

    call monolis_check_inside_in_aabb(coord, BB, ths, is_inside)

    call monolis_test_check_eq_L1("monolis_check_inside_in_aabb 3", is_inside, .false.)


    coord(1) = 2.0d0
    coord(2) = 20.0d0
    coord(3) = 0.0d0

    call monolis_check_inside_in_aabb(coord, BB, ths, is_inside)

    call monolis_test_check_eq_L1("monolis_check_inside_in_aabb 4", is_inside, .false.)
  end subroutine monolis_check_inside_in_aabb_test

  !> main test subroutine
  subroutine monolis_utils_kdtree_test()
    implicit none

    call monolis_kdtree_init_test()
    call monolis_kdtree_overlap_test()

    call monolis_std_global_log_string("monolis_kdtree_finalize_main")
    call monolis_std_global_log_string("monolis_kdtree_get_BB_including_coordinates_main")
    call monolis_std_global_log_string("monolis_kdtree_init_by_BB_main")
  end subroutine monolis_utils_kdtree_test

  !> unit test
  subroutine monolis_kdtree_init_test()
    implicit none
    type(monolis_kdtree_structure) :: monolis_kdtree
    integer(kint) :: n_BB
    integer(kint) :: BB_id(5)
    real(kdouble) :: BB(6,5), pos(3)
    integer(kint), allocatable :: ids(:)

    call monolis_std_global_log_string("monolis_kdtree_init_by_BB")
    call monolis_std_global_log_string("monolis_kdtree_get_BB_including_coordinates")
    call monolis_std_global_log_string("monolis_kdtree_finalize")

    !> initialize
    BB_id(1) = 10
    BB_id(2) = 20
    BB_id(3) = 30
    BB_id(4) = 40
    BB_id(5) = 50

    BB(1,1) =  1.0d0; BB(2,1) =  2.0d0; BB(3,1) = 51.0d0; BB(4,1) = 52.0d0; BB(5,1) = 101.0d0; BB(6,1) = 102.0d0;
    BB(1,2) = 11.0d0; BB(2,2) = 12.0d0; BB(3,2) = 61.0d0; BB(4,2) = 62.0d0; BB(5,2) = 111.0d0; BB(6,2) = 112.0d0;
    BB(1,3) = 21.0d0; BB(2,3) = 22.0d0; BB(3,3) = 71.0d0; BB(4,3) = 72.0d0; BB(5,3) = 121.0d0; BB(6,3) = 122.0d0;
    BB(1,4) = 31.0d0; BB(2,4) = 32.0d0; BB(3,4) = 81.0d0; BB(4,4) = 82.0d0; BB(5,4) = 131.0d0; BB(6,4) = 132.0d0;
    BB(1,5) = 41.0d0; BB(2,5) = 42.0d0; BB(3,5) = 91.0d0; BB(4,5) = 92.0d0; BB(5,5) = 141.0d0; BB(6,5) = 142.0d0;

    n_BB = 5

    call monolis_kdtree_init_by_BB(monolis_kdtree, n_BB, BB_id, BB)

    !> case 1
    pos(1) = 1.0d0
    pos(2) = 1.0d0
    pos(3) = 1.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 1", n_BB, 0)

    !> case 2
    pos(1) = 1.0d0
    pos(2) = 51.0d0
    pos(3) = 101.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 2", n_BB, 1)
    call monolis_test_check_eq_I1("monolis_kdtree_init_test 2", ids(1), 10)

    !> case 3
    pos(1) = 11.0d0
    pos(2) = 61.0d0
    pos(3) = 111.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 3", n_BB, 1)
    call monolis_test_check_eq_I1("monolis_kdtree_init_test 3", ids(1), 20)

    !> case 4
    pos(1) = 21.0d0
    pos(2) = 71.0d0
    pos(3) = 121.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 4", n_BB, 1)
    call monolis_test_check_eq_I1("monolis_kdtree_init_test 4", ids(1), 30)

    !> case 5
    pos(1) = 31.0d0
    pos(2) = 81.0d0
    pos(3) = 131.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 5", n_BB, 1)
    call monolis_test_check_eq_I1("monolis_kdtree_init_test 5", ids(1), 40)

    !> case 6
    pos(1) = 41.0d0
    pos(2) = 91.0d0
    pos(3) = 141.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)

    call monolis_test_check_eq_I1("monolis_kdtree_init_test 6", n_BB, 1)
    call monolis_test_check_eq_I1("monolis_kdtree_init_test 6", ids(1), 50)

    call monolis_kdtree_finalize(monolis_kdtree)
  end subroutine monolis_kdtree_init_test

  subroutine monolis_kdtree_overlap_test()
    implicit none
    type(monolis_kdtree_structure) :: monolis_kdtree
    integer(kint) :: n_BB
    integer(kint) :: BB_id(5), ids_ans(5)
    real(kdouble) :: BB(6,5), pos(3)
    integer(kint), allocatable :: ids(:)

    call monolis_std_global_log_string("monolis_kdtree_init_by_BB")
    call monolis_std_global_log_string("monolis_kdtree_get_BB_including_coordinates")
    call monolis_std_global_log_string("monolis_kdtree_finalize")

    !> initialize
    BB_id(1) = 50
    BB_id(2) = 40
    BB_id(3) = 30
    BB_id(4) = 20
    BB_id(5) = 10

    BB(1,1) = 0.0d0; BB(2,1) = 5.0d0; BB(3,1) = 0.0d0; BB(4,1) = 5.0d0; BB(5,1) = 0.0d0; BB(6,1) = 5.0d0;
    BB(1,2) = 0.0d0; BB(2,2) = 4.0d0; BB(3,2) = 0.0d0; BB(4,2) = 4.0d0; BB(5,2) = 0.0d0; BB(6,2) = 4.0d0;
    BB(1,3) = 0.0d0; BB(2,3) = 3.0d0; BB(3,3) = 0.0d0; BB(4,3) = 3.0d0; BB(5,3) = 0.0d0; BB(6,3) = 3.0d0;
    BB(1,4) = 0.0d0; BB(2,4) = 2.0d0; BB(3,4) = 0.0d0; BB(4,4) = 2.0d0; BB(5,4) = 0.0d0; BB(6,4) = 2.0d0;
    BB(1,5) = 0.0d0; BB(2,5) = 1.0d0; BB(3,5) = 0.0d0; BB(4,5) = 1.0d0; BB(5,5) = 0.0d0; BB(6,5) = 1.0d0;

    n_BB = 5

    call monolis_kdtree_init_by_BB(monolis_kdtree, n_BB, BB_id, BB)

    !> case 1
    pos(1) = 1.0d0
    pos(2) = 1.0d0
    pos(3) = 1.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)
    call monolis_qsort_I_1d(ids, 1, 5)

    ids_ans(1) = 10
    ids_ans(2) = 20
    ids_ans(3) = 30
    ids_ans(4) = 40
    ids_ans(5) = 50

    call monolis_test_check_eq_I1("monolis_kdtree_overlap_test 1", n_BB, 5)
    call monolis_test_check_eq_I("monolis_kdtree_overlap_test 1", ids, ids_ans)

    !> case 2
    pos(1) = 2.0d0
    pos(2) = 2.0d0
    pos(3) = 2.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)
    call monolis_qsort_I_1d(ids, 1, 4)

    ids_ans(1) = 20
    ids_ans(2) = 30
    ids_ans(3) = 40
    ids_ans(4) = 50

    call monolis_test_check_eq_I1("monolis_kdtree_overlap_test 2", n_BB, 4)
    call monolis_test_check_eq_I("monolis_kdtree_overlap_test 2", ids, ids_ans(1:4))

    !> case 3
    pos(1) = 3.0d0
    pos(2) = 3.0d0
    pos(3) = 3.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)
    call monolis_qsort_I_1d(ids, 1, 3)

    ids_ans(1) = 30
    ids_ans(2) = 40
    ids_ans(3) = 50

    call monolis_test_check_eq_I1("monolis_kdtree_overlap_test 3", n_BB, 3)
    call monolis_test_check_eq_I("monolis_kdtree_overlap_test 3", ids, ids_ans(1:3))

    !> case 4
    pos(1) = 4.0d0
    pos(2) = 4.0d0
    pos(3) = 4.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)
    call monolis_qsort_I_1d(ids, 1, 2)

    ids_ans(1) = 40
    ids_ans(2) = 50

    call monolis_test_check_eq_I1("monolis_kdtree_overlap_test 4", n_BB, 2)
    call monolis_test_check_eq_I("monolis_kdtree_overlap_test 4", ids, ids_ans(1:2))

    !> case 5
    pos(1) = 5.0d0
    pos(2) = 5.0d0
    pos(3) = 5.0d0

    call monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, ids)
    call monolis_qsort_I_1d(ids, 1, 1)

    ids_ans(1) = 50

    call monolis_test_check_eq_I1("monolis_kdtree_overlap_test 5", n_BB, 1)
    call monolis_test_check_eq_I("monolis_kdtree_overlap_test 5", ids, ids_ans(1:1))

    call monolis_kdtree_finalize(monolis_kdtree)
  end subroutine monolis_kdtree_overlap_test
end module mod_gg_tools_def_bucket_test
