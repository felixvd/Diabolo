import numpy as np

"""
This file contains a environment for the 2D diabolo simulator:
The equations of motions of the diabolo can be viewed on 
https://www.overleaf.com/project/60c7b833bac7280c1b660be8
"""


def update_position(pos_t_d, dot_pos_t_d, pos_t_s, l, wt, deltat=0.00001, g=9.8):
    """
    pos_t_d: a tuple with two elements, (diabolo_pos_xt, diabolo_pos_yt)
    dot_pos_t_d: a tuple with two elements, (dot_diabolo_pos_xt, dot_diabolo_pos_yt)
    pos_t_s: a tuple with four elements, (left_stick_pos_xt, left_stick_pos_yt, right_stick_pos_xt, right_stick_pos_yt)
    """
    xt_d, yt_d = pos_t_d
    dot_xt_d, dot_yt_d = dot_pos_t_d
    xt_l_s, yt_l_s, xt_r_s, yt_r_s = pos_t_s 

    # Unconstrained position, with velocity-verlet euler integration
    xtplus1_d = xt_d + dot_xt_d * deltat
    ytplus1_d = yt_d + dot_yt_d * deltat - 0.5 * g * deltat**2

    # Constrain Position due to string
    ## compute the center of the two sticks at time t
    x_center = (xt_l_s + xt_r_s) / 2.
    y_center = (yt_l_s + yt_r_s) / 2.
    ## Compute the parameters of the line
    m = (ytplus1_d-y_center)/(xtplus1_d-x_center)
    n = -1
    p = y_center - m * x_center
    ## Compute the parameters of the ellipse define by the two stick tips
    a = l/2.
    b = np.sqrt(a**2-(np.abs(xt_l_s-xt_r_s)/2)**2)
    
    if a**2*m**2+b**2*n**2 != 0 and n!=0 and a*b != 0:




