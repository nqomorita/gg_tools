program gg_tools_test
  use mod_monolis_utils
  use mod_gg_tools
  use mod_gg_tools_def_bucket_test
  use mod_gg_tools_distance_determination_test
  implicit none

  call monolis_mpi_initialize()

  call monolis_utils_aabb_test()
  call monolis_utils_kdtree_test()
  call gg_tools_def_bucket_test()
  call gg_tools_distance_determination_test()

  call monolis_mpi_finalize()
end program gg_tools_test
