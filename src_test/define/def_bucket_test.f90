!> バケットセルテストモジュール
module mod_gg_tools_def_bucket_test
  use mod_gg_tools
  use mod_monolis_utils
  implicit none

contains

  subroutine gg_tools_def_bucket_test()
    implicit none
  end subroutine gg_tools_def_bucket_test

  subroutine gg_tools_dlb_initialize_test()
    implicit none
    call monolis_std_log_string("gg_tools_dlb_initialize")
  end subroutine gg_tools_dlb_initialize_test

  subroutine gg_tools_dlb_finalize_test()
    implicit none
    call monolis_std_log_string("gg_tools_dlb_finalize")
  end subroutine gg_tools_dlb_finalize_test
end module mod_gg_tools_def_bucket_test
