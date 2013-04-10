  module typedefs
    implicit none
  contains
    type prim_update_options
       integer :: algorithm_choice
       type(prim_update_FV_options) :: FV_opts
    end type prim_update_options
    type prim_update_FV_options
       integer :: ghost_pts
       integer :: spatial_order
       integer :: grid_motion
    end type prim_update_FV_options
  end module typedefs
