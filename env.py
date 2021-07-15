import numpy as np

"""
This file contains a environment for the 2D diabolo simulator:
The equations of motions of the diabolo can be viewed on 
https://www.overleaf.com/project/60c7b833bac7280c1b660be8
"""


def update_position(pos_t_d, dot_pos_t_d, pos_t_s, dot_pos_t_s, l, wt, deltat=0.00001, g=9.8):
    """
    Update the position of the diaoblo and the stick tips

    Args:
        pos_t_d: the position of the diabolo at time t
            a tuple with two elements, (diabolo_pos_xt, diabolo_pos_yt)
        dot_pos_t_d: the velocity of the diabolo at time t, 
            a tuple with two elements, (dot_diabolo_pos_xt, dot_diabolo_pos_yt)
        pos_t_s: the positions of the stick tip at time t
            a tuple with four elements, (left_stick_pos_xt, left_stick_pos_yt, right_stick_pos_xt, right_stick_pos_yt)
        dot_pos_t_s: the velocity of the stick tip at time t
            a tuple with four elements, (dot_left_stick_pos_xt, dot_left_stick_pos_yt, dot_right_stick_pos_xt, dot_right_stick_pos_yt)
        l: the length of the string, a scalar
        wt: the rotational speed of the diabolo, a scalar
        deltat: the time interval for the velocity verlet algorithm
        g: the gravitional constant

    Return:
        pos_t_d
        dot_pos_t_d
        pos_t_s
        dot_pos_t_s
        wt
    """
    xt_d, yt_d = pos_t_d
    dot_xt_d, dot_yt_d = dot_pos_t_d
    xt_l_s, yt_l_s, xt_r_s, yt_r_s = pos_t_s 

    # Step1: Unconstrained position, with velocity-verlet euler integration
    xtplus1_d = xt_d + dot_xt_d * deltat
    ytplus1_d = yt_d + dot_yt_d * deltat - 0.5 * g * deltat**2

    # Step2: Constrain Position due to string
    ## Compute the center of the two sticks at time t
    x_center = (xt_l_s + xt_r_s) / 2.
    y_center = (yt_l_s + yt_r_s) / 2.
    ## Compute the parameters of the line
    m = (ytplus1_d-y_center)/(xtplus1_d-x_center)
    n = -1
    p = y_center - m * x_center
    ## Compute the parameters of the ellipse define by the two stick tips
    a = l/2.
    b = np.sqrt(a**2-(np.abs(xt_l_s-xt_r_s)/2)**2)
    ## Project the new position to the tangent plane of the ellipse
    if a**2*m**2+b**2*n**2 != 0 and n!=0 and a*b != 0:
        x_tmp1 = np.sqrt(a**2*b**2*n**2*(a**2*m**2+b**2*n**2-p**2))+a**2*m*p
        x_tmp1 = x_tmp1 / (a**2*m**2+b**2*n**2)
        y_tmp1 = m*np.sqrt(a**2*b**2*n**2*(a**2*m**2+b**2*n**2-p**2))-b**2*n**2*p
        y_tmp1 = y_tmp1 / (a**2*m**2*n+b**2*n**3)

        x_tmp2 = np.sqrt(a**2*b**2*n**2*(a**2*m**2+b**2*n**2-p**2))-a**2*m*p
        x_tmp2 = x_tmp2 / (a**2*m**2+b**2*n**2)
        y_tmp2 = -(m*np.sqrt(a**2*b**2*n**2*(a**2*m**2+b**2*n**2-p**2))+b**2*n**2*p)
        y_tmp2 = y_tmp2 / (a**2*m**2*n+b**2*n**3)

        if x_tmp1*xtplus1_d > 0 and y_tmp1*ytplus1_d > 0:
            x_c, y_c = x_tmp1, y_tmp1
        else:
            x_c, y_c = x_tmp2, y_tmp2

    if n == 0 and m != 0 and a != 0 and b!=0:
        x_c = p / m
        y_tmp1 = b*np.sqrt(a**2-p**2/m**2)/a
        y_tmp2 = -b*np.sqrt(a**2-p**2/m**2)/a
        
        if y_tmp1*ytplus1_d > 0:
            y_c = y_tmp1
        else:
            y_c = y_tmp2

    # Step3: Add pull velocity if the sticks are moving
    if np.linalg.norm(dot_pos_t_s) > 1e-2:
        v_pull = (x_c - xtplus1_d, y_c - ytplus1_d)
        dot_pos_t_d = (dot_pos_t_d[0]+v_pull[0], dot_pos_t_d[1]+v_pull[1])

    # Step4: Constrain velocity
    u = -b**2*c_c / (a**2*y_c)
    v = np.array(dot_pos_t_d)
    n = np.array([0, y_c-u*x_c])
    if np.cross(v, n) < 0:
        cos_vn = np.inner(v,n)/(np.linalg.norm(v)*np.linalg.norm(n))
        v_c = v * cos_vn
        dot_pos_t_d = v_c
    
    # Update rotational velocity 

    return pos_t_d, dot_pos_t_d, pos_t_s, dot_pos_t_s, wt 
