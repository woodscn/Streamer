import numpy as np
def make_options(grid_motion,spatial_order):
    nopts = 1000
    out = np.zeros(nopts)
    grid_motions = {'eulerian':0,
                    'lagrangian':1,
                    'half-n-half':2,
                    'grid_angle_preserving':3}
    spatial_orders = {1:1,2:2}
    out[0] = 1
    out[100] = 1
    out[101] = spatial_orders[spatial_order]
    out[102] = grid_motions[grid_motion]
    out[103] = 0
    return out
    
