#Options meanings
# [1-2]: controls which prim_update algorithm to use
# [3-5]: sets values for dxi, deta, dzeta (0=>1., 1=>0.5, 2=>0.25)
# [6-7]: controls grid motion. 
# 6: 0=>Eulerian, 1=>Lagrangian-esque, 2=>Constant, 3=>Angle-preserving
# 7: h0 = (1=>.25, 2=>.5, 3=>.999)
# [101]: reports how many boundary ghost points are present
# [102]: controls spatial order of accuracy
# [104]: Controls type of time step (constant or CFL)
# [301]: Controls type of grid motion
