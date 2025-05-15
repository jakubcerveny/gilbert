#!/usr/bin/env python3
#
# To the extent possible under law, the person who associated CC0 with
# this project has waived all copyright and related or neighboring rights
# to this project.
# 
# You should have received a copy of the CC0 legalcode along with this
# work. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
#

# This is a accompanying reference implementation for the
# Gilbert2D and Gilbert3D functions referenced in the
# paper https://github.com/jakubcerveny/gilbert-paper
#
# The code here's main focus is legibility, not execution
# speed or memory usage.
# There are both asynchronous functions that closely match
# the psuedo code from the paper as well as synchronous
# functions with their corresponding d2xyz and xyz2d
# functions.
#
# Running from the command line will list options.
#

VERBOSE = 1

ADAPT_METHOD = {
    "HARMONY" : 0,
    "HAMILTONIAN": 1,
    "AXIS": 2
}

###################################################################
#    __       __               ___              __  _             
#   / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
#  / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
# /_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
#           /_/                                                   
###################################################################


def sgn(t): return 1 if (t>0) else (-1 if (t<0) else 0)

# vector element wise absolute value
#
def v_delta(v): return list(map( lambda t: sgn(t), v ))
def v_divq(v,q): return list(map( lambda t: int(t/q), v ))
def v_div2(v): return v_divq(v,2)

def v_add(u,v):
    w = []
    for i in range(len(u)): w.append( u[i] + v[i] )
    return w

def v_neg(u):
    w = []
    for i in range(len(u)): w.append( -u[i] )
    return w

def v_mul(c,v):
    w = []
    for i in range(len(v)): w.append( c * v[i] )
    return w

def dot_v(u,v):
    s = 0
    for i in range(len(u)): s += u[i] * v[i]
    return s

def abs_sum_v(v): return sum(list(map(lambda t: abs(t), v)))


# divide by two, force even
# adding one if necessary
#
def d2e(t):
    if t == 0: return 0

    T = abs(t)
    s = sgn(t)
    T2 = int(T/2)
    if (T2%2) == 0: return s*T2
    return s*(T2 + 1)

# divide by two, force odd
# adding one if necessary
#
def d2u(t):
    if t == 0: return 0

    T = abs(t)
    s = sgn(t)
    T2 = int(T/2)
    if (T2%2) == 1: return s*T2
    return s*(T2 + 1)


# divide by q, force even
# adding one if necessary
#
def dqe(t,q):
    if t == 0: return 0

    T = abs(t)
    s = sgn(t)
    Tq = int(T/q)
    if (Tq%2) == 0: return s*Tq
    return s*(Tq + 1)

# divide by q, force odd
# adding one if necessary
#
def dqu(t,q):
    if t == 0: return 0

    T = abs(t)
    s = sgn(t)
    Tq = int(T/q)
    if (Tq%2) == 1: return s*Tq
    return s*(Tq + 1)

def inBounds(q, p, a, b, g=None):
    _a = [a[0],a[1],0]
    _b = [b[0],b[1],0]
    _g = [0,0,0]

    _p = [p[0],p[1],0]
    _q = [q[0],q[1],0]

    if len(p) > 2: _p[2] = p[2]
    if len(q) > 2: _q[2] = q[2]

    if len(a) > 2: _a[2] = a[2]
    if len(b) > 2: _b[2] = b[2]

    if g is None:
        if   (_a[0] == 0) and (_b[0] == 0): _g[0] = 1
        elif (_a[1] == 0) and (_b[1] == 0): _g[1] = 1
        elif (_a[2] == 0) and (_b[2] == 0): _g[2] = 1
    else:
        _g = [g[0], g[1], g[2]]

    for xyz in range(3):
        d = _a[xyz] + _b[xyz] + _g[xyz]
        if d < 0:
            if ((_q[xyz] > _p[xyz]) or
                (_q[xyz] <= (_p[xyz] + d))): return False
        else:
            if ((_q[xyz] < _p[xyz]) or
                (_q[xyz] >= (_p[xyz] + d))): return False
    return True

def Hilbert2x2x2Async(p, alpha, beta, gamma):
    u = [p[0], p[1], p[2]]

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)
    d_gamma = v_delta(gamma)

    yield [ u[0], u[1], u[2] ]

    u = v_add(u, d_beta)
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, d_gamma)
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, v_neg(d_beta))
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, d_alpha)
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, d_beta)
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, v_neg(d_gamma))
    yield [ u[0], u[1], u[2] ]

    u = v_add(u, v_neg(d_beta))
    yield [ u[0], u[1], u[2] ]

def Hilbert2x2x2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma):
    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)
    d_gamma = v_delta(gamma)

    d_idx = dst_idx - cur_idx

    if    d_idx == 0: return [p[0], p[1], p[2]]
    elif  d_idx == 1: return v_add(p, d_beta)
    elif  d_idx == 2: return v_add(p, v_add(d_beta, d_gamma))
    elif  d_idx == 3: return v_add(p, d_gamma)
    elif  d_idx == 4: return v_add(p, v_add(d_alpha, d_gamma))
    elif  d_idx == 5: return v_add(p, v_add(d_alpha, v_add(d_beta, d_gamma)))
    elif  d_idx == 6: return v_add(p, v_add(d_alpha, d_beta))
    elif  d_idx == 7: return v_add(p, d_alpha)

    return [-1,-1,-1]


###################################################################
#    __       __               ___              __  _             
#   / /  ___ / /__  ___ ____  / _/_ _____  ____/ /_(_)__  ___  ___
#  / _ \/ -_) / _ \/ -_) __/ / _/ // / _ \/ __/ __/ / _ \/ _ \(_-<
# /_//_/\__/_/ .__/\__/_/   /_/ \_,_/_//_/\__/\__/_/\___/_//_/___/
#           /_/                                                   
###################################################################


####################################################
#   ______ ____           __  ____ ___    __    __ 
#  / ___(_) / /  ___ ____/ /_|_  // _ \__/ /___/ /_
# / (_ / / / _ \/ -_) __/ __//_ </ // /_  __/_  __/
# \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
#                                                  
####################################################

########################
#  __ _ ____  _ _ _  __ 
# / _` (_-< || | ' \/ _|
# \__,_/__/\_, |_||_\__|
#          |__/         
########################


def Gilbert3DS0Async(p, alpha, beta, gamma):
    alpha2  = v_div2(alpha)
    d_alpha = v_delta(alpha)

    a   = abs_sum_v(alpha)
    a2  = abs_sum_v(alpha2)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)

    yield from Gilbert3DAsync( p, alpha2, beta, gamma )
    yield from Gilbert3DAsync( v_add(p, alpha2),
                               v_add(alpha, v_neg(alpha2)), beta, gamma )

def Gilbert3DS1Async(p, alpha, beta, gamma):
    alpha2 = v_div2(alpha)
    gamma3 = v_divq(gamma, 3)

    d_alpha = v_delta(alpha)
    d_gamma = v_delta(gamma)

    a = abs_sum_v(alpha)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    g3 = abs_sum_v(gamma3)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)

    if (g > 2) and ((g3 % 2) == 1):
        gamma3 = v_add(gamma3, d_gamma)

    yield from Gilbert3DAsync( p,
                               gamma3,
                               alpha2,
                               beta )

    yield from Gilbert3DAsync( v_add(p, gamma3),
                               alpha,
                               beta,
                               v_add(gamma, v_neg(gamma3)) )

    yield from Gilbert3DAsync( v_add(p,
                                     v_add( v_add(alpha , v_neg(d_alpha)),
                                            v_add(gamma3, v_neg(d_gamma)) ) ),
                               v_neg(gamma3),
                               v_neg( v_add(alpha, v_neg(alpha2)) ),
                               beta )


def Gilbert3DS2Async(p, alpha, beta, gamma):
    alpha2 = v_div2(alpha)
    beta3  = v_divq(beta, 3)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)

    a2 = abs_sum_v(alpha2)
    b3 = abs_sum_v(beta3)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)

    if (b > 2) and ((b3 % 2) == 1):
        beta3 = v_add(beta3, d_beta)

    yield from Gilbert3DAsync( p,
                               beta3,
                               gamma,
                               alpha2 )

    yield from Gilbert3DAsync( v_add(p, beta3),
                               alpha,
                               v_add(beta, v_neg(beta3)),
                               gamma )

    yield from Gilbert3DAsync( v_add(p,
                                     v_add( v_add(alpha, v_neg(d_alpha)),
                                            v_add(beta3, v_neg(d_beta)) ) ),
                               v_neg(beta3),
                               gamma,
                               v_neg( v_add(alpha, v_neg(alpha2)) ) )



def Gilbert3DJ0Async(p, alpha, beta, gamma):
    alpha2 = v_div2(alpha)
    beta2  = v_div2(beta)
    gamma2 = v_div2(gamma)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)
    d_gamma = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 1): alpha2 = v_add(alpha2, d_alpha)
    if (b > 2) and ((b2 % 2) == 1): beta2  = v_add(beta2, d_beta)
    if (g > 2) and ((g2 % 2) == 1): gamma2 = v_add(gamma2, d_gamma)

    yield from Gilbert3DAsync( p,
                               beta2, gamma2, alpha2 )

    yield from Gilbert3DAsync( v_add(p, beta2),
                               gamma, alpha2, v_add(beta, v_neg(beta2)) )

    yield from Gilbert3DAsync( v_add( p, v_add( v_add(beta2, v_neg(d_beta)), v_add(gamma, v_neg(d_gamma)) ) ),
                               alpha, v_neg(beta2), v_neg( v_add(gamma, v_neg(gamma2)) ) )


    yield from Gilbert3DAsync( v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_add(gamma, v_neg(d_gamma))) ) ),
                               v_neg(gamma), v_neg( v_add(alpha, v_neg(alpha2)) ), v_add(beta, v_neg(beta2)) )

    yield from Gilbert3DAsync( v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_neg(d_beta)) ) ),
                               v_neg(beta2), gamma2, v_neg( v_add(alpha, v_neg(alpha2)) ) )


def Gilbert3DJ1Async(p, alpha, beta, gamma):
    alpha2 = v_div2(alpha)
    beta2  = v_div2(beta)
    gamma2 = v_div2(gamma)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)
    d_gamma = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 0): alpha2 = v_add(alpha2, d_alpha)
    if (b > 2) and ((b2 % 2) == 1): beta2  = v_add(beta2, d_beta)
    if (g > 2) and ((g2 % 2) == 1): gamma2 = v_add(gamma2, d_gamma)

    yield from Gilbert3DAsync( p,
                               gamma2, alpha2, beta2 )

    yield from Gilbert3DAsync( v_add(p, gamma2),
                               beta, v_add(gamma, v_neg(gamma2)), alpha2 )


    yield from Gilbert3DAsync( v_add( p, v_add( v_add(gamma2, v_neg(d_gamma)), v_add(beta, v_neg(d_beta)) ) ),
                               alpha, v_neg(v_add(beta, v_neg(beta2))), v_neg(gamma2) )


    yield from Gilbert3DAsync( v_add( p , v_add( v_add(alpha, v_neg(d_alpha)), v_add( v_add(beta, v_neg(d_beta)), gamma2 ) ) ),
                               v_neg(beta), v_add(gamma, v_neg(gamma2)), v_neg(v_add(alpha, v_neg(alpha2))) )


    yield from Gilbert3DAsync( v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(gamma2, v_neg(d_gamma)) ) ),
                               v_neg(gamma2), v_neg(v_add(alpha, v_neg(alpha2))), beta2 )


def Gilbert3DJ2Async(p, alpha, beta, gamma):
    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)
    gamma2  = v_div2(gamma)

    d_alpha  = v_delta(alpha)
    d_beta   = v_delta(beta)
    d_gamma  = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 0): alpha2 = v_add(alpha2, d_alpha)
    if (b > 2) and ((b2 % 2) == 1): beta2  = v_add(beta2, d_beta)
    if (g > 2) and ((g2 % 2) == 1): gamma2 = v_add(gamma2, d_gamma)

    yield from Gilbert3DAsync( p,
                               beta2, gamma, alpha2 )

    yield from Gilbert3DAsync( v_add(p, beta2),
                               gamma2, alpha, v_add(beta, v_neg(beta2)) )

    yield from Gilbert3DAsync( v_add(p, v_add(beta2, gamma2)),
                               alpha, v_add(beta, v_neg(beta2)), v_add(gamma, v_neg(gamma2)) )

    yield from Gilbert3DAsync( v_add( p, v_add( v_add( alpha, v_neg(d_alpha) ), v_add( v_add(beta2, v_neg(d_beta)), gamma2 ) ) ),
                               v_neg(beta2), v_add(gamma, v_neg(gamma2)), v_neg(v_add(alpha, v_neg(alpha2))) )

    yield from Gilbert3DAsync( v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add( gamma2, v_neg(d_gamma) ) ) ),
                               v_neg(gamma2), v_neg(v_add(alpha, v_neg(alpha2))), beta2)


def Gilbert3DAsync(p, alpha, beta, gamma):
    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    (a0,b0,g0) = ( a % 2, b % 2, g % 2 )

    if VERBOSE > 0: print("#G3Dasync: p:", p, "alpha:", alpha, "beta:", beta, "gamma:", gamma, "(", a0,b0,g0,")")

    if (a == 2) and (b == 2) and (g == 2):

        if VERBOSE > 0: print("#  G3D.h2x2x2")

        yield from Hilbert2x2x2Async(p, alpha, beta, gamma)
        return

    if a == 1:

        if VERBOSE > 0: print("#  G2D.bg")

        yield from Gilbert2DAsync(p, beta, gamma)
        return

    if b == 1:

        if VERBOSE > 0: print("#  G2D.ag")

        yield from Gilbert2DAsync(p, alpha, gamma)
        return

    if g == 1:

        if VERBOSE > 0: print("#  G2D.ab")

        yield from Gilbert2DAsync(p, alpha, beta)
        return

    if (((3*a) > (5*b)) and
        ((3*a) > (5*g))):

        if VERBOSE > 0: print("#  G3D.S0")

        yield from Gilbert3DS0Async(p, alpha, beta, gamma)
        return

    if (((2*b) > (3*g)) or
        ((2*b) > (3*a))):

        if VERBOSE > 0: print("#  G3D.S2")

        yield from Gilbert3DS2Async(p, alpha, beta, gamma)
        return

    if ((2*g) > (3*b)):

        if VERBOSE > 0: print("#  G3D.S1")

        yield from Gilbert3DS1Async(p, alpha, beta, gamma)
        return

    if g0 == 0:

        if VERBOSE > 0: print("#  G3D.J0")

        yield from Gilbert3DJ0Async(p, alpha, beta, gamma)
        return

    if (a0 == 0) or (b0 == 0):

        if VERBOSE > 0: print("#  G3D.J1")

        yield from Gilbert3DJ1Async(p, alpha, beta, gamma)
        return

    if VERBOSE > 0: print("#  G3D.J2")

    yield from Gilbert3DJ2Async(p, alpha, beta, gamma)

########################
#  __ _ ____  _ _ _  __ 
# / _` (_-< || | ' \/ _|
# \__,_/__/\_, |_||_\__|
#          |__/         
########################


########################
#              ___    _ 
# __ ___  _ __|_  )__| |
# \ \ / || |_ // // _` |
# /_\_\\_, /__/___\__,_|
#      |__/             
########################


########################
#              ___    _ 
# __ ___  _ __|_  )__| |
# \ \ / || |_ // // _` |
# /_\_\\_, /__/___\__,_|
#      |__/             
########################


########################
#     _ ___             
#  __| |_  )_ ___  _ ___
# / _` |/ /\ \ / || |_ /
# \__,_/___/_\_\\_, /__|
#               |__/    
########################

def Gilbert3DS0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma):

    alpha2  = v_div2(alpha)
    d_alpha = v_delta(alpha)

    a   = abs_sum_v(alpha)
    a2  = abs_sum_v(alpha2)

    b   = abs_sum_v(beta)
    g   = abs_sum_v(gamma)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)

    nxt_idx = cur_idx + (a2*b*g)
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                alpha2, beta, gamma )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add(p, alpha2),
                            v_add(alpha, v_neg(alpha2)), beta, gamma )


def Gilbert3DS1_d2xyz( dst_idx, cur_idx, p, alpha, beta, gamma ):
    alpha2 = v_div2(alpha)
    gamma3 = v_divq(gamma, 3)

    d_alpha = v_delta(alpha)
    d_gamma = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    g3 = abs_sum_v(gamma3)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)

    if (g > 2) and ((g3 % 2) == 1):
        gamma3 = v_add(gamma3, d_gamma)
        g3 = abs_sum_v(gamma3)

    nxt_idx = cur_idx + (g3*a2*b)
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                gamma3, alpha2, beta )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*b*(g-g3))
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add(p, gamma3),
                                alpha, beta, v_add(gamma, v_neg(gamma3)) )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add(p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(gamma3, v_neg(d_gamma)) ) ),
                            v_neg(gamma3), v_neg(v_add(alpha, v_neg(alpha2))), beta )

def Gilbert3DS2_d2xyz( dst_idx, cur_idx, p, alpha, beta, gamma ):
    alpha2  = v_div2(alpha)
    beta3   = v_divq(beta, 3)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b3 = abs_sum_v(beta3)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)

    if (b > 2) and((b3 % 2) == 1):
        beta3 = v_add(beta3, d_beta)
        b3 = abs_sum_v(beta3)


    nxt_idx = cur_idx + (b3*g*a2)
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                beta3, gamma, alpha2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*(b-b3)*g)
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add(p, beta3),
                                alpha, v_add(beta, v_neg(beta3)), gamma )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta3, v_neg(d_beta)) ) ),
                            v_neg(beta3), gamma, v_neg(v_add(alpha, v_neg(alpha2))) )


def Gilbert3DJ0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma):
    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)
    gamma2  = v_div2(gamma)

    d_alpha  = v_delta(alpha)
    d_beta   = v_delta(beta)
    d_gamma  = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 1):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)
    if (b > 2) and ((b2 % 2) == 1):
        beta2  = v_add(beta2, d_beta)
        b2 = abs_sum_v(beta2)
    if (g > 2) and ((g2 % 2) == 1):
        gamma2 = v_add(gamma2, d_gamma)
        g2 = abs_sum_v(gamma2)

    nxt_idx = cur_idx + (b2*g2*a2)
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                beta2, gamma2, alpha2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (g*a2*(b-b2))
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add(p, beta2),
                                gamma, alpha2, v_add(beta, v_neg(beta2)) )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*b2*(g-g2))
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p, v_add( v_add(beta2, v_neg(d_beta)), v_add(gamma, v_neg(d_gamma)) ) ),
                                alpha, v_neg(beta2), v_neg( v_add(gamma, v_neg(gamma2)) ) )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (g*(a-a2)*(b-b2))
    if ((cur_idx <= dst_idx) and (dst_idx < nxt_idx)):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_add(gamma, v_neg(d_gamma))) ) ),
                                v_neg(gamma), v_neg( v_add(alpha, v_neg(alpha2)) ), v_add(beta, v_neg(beta2)) )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_neg(d_beta)) ) ),
                            v_neg(beta2), gamma2, v_neg( v_add(alpha, v_neg(alpha2)) ) )

def Gilbert3DJ1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma):

    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)
    gamma2  = v_div2(gamma)

    d_alpha  = v_delta(alpha)
    d_beta   = v_delta(beta)
    d_gamma  = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 0):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)
    if (b > 2) and ((b2 % 2) == 1):
        beta2  = v_add(beta2, d_beta)
        b2 = abs_sum_v(beta2)
    if (g > 2) and ((g2 % 2) == 1):
        gamma2 = v_add(gamma2, d_gamma)
        g2 = abs_sum_v(gamma2)

    nxt_idx = cur_idx + (g2*a2*b2)
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                gamma2, alpha2, beta2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (b*(g-g2)*a2)
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p, gamma2 ),
                                beta, v_add(gamma, v_neg(gamma2)), alpha2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*(b-b2)*g2)
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p, v_add( v_add(gamma2, v_neg(d_gamma)), v_add(beta, v_neg(d_beta)) ) ),
                                alpha, v_neg(v_add(beta, v_neg(beta2))), v_neg(gamma2) )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (b*(g-g2)*(a-a2))
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p , v_add( v_add(alpha, v_neg(d_alpha)), v_add( v_add(beta, v_neg(d_beta)), gamma2 ) ) ),
                                v_neg(beta), v_add(gamma, v_neg(gamma2)), v_neg(v_add(alpha, v_neg(alpha2))) )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(gamma2, v_neg(d_gamma)) ) ),
                            v_neg(gamma2), v_neg(v_add(alpha, v_neg(alpha2))), beta2 )


def Gilbert3DJ2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma):

    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)
    gamma2  = v_div2(gamma)

    d_alpha  = v_delta(alpha)
    d_beta   = v_delta(beta)
    d_gamma  = v_delta(gamma)

    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)
    g2 = abs_sum_v(gamma2)

    if (a > 2) and ((a2 % 2) == 0):
        alpha2 = v_add(alpha2, d_alpha)
        a2 = abs_sum_v(alpha2)
    if (b > 2) and ((b2 % 2) == 1):
        beta2  = v_add(beta2, d_beta)
        b2 = abs_sum_v(beta2)
    if (g > 2) and ((g2 % 2) == 1):
        gamma2 = v_add(gamma2, d_gamma)
        g2 = abs_sum_v(gamma2)

    nxt_idx = cur_idx + (b2*g*a2)
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                p,
                                beta2, gamma, alpha2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (g2*a*(b-b2))
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add(p, beta2),
                                gamma2, alpha, v_add(beta, v_neg(beta2)) )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*(b-b2)*(g-g2))
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add(p, v_add(beta2, gamma2)),
                                alpha, v_add(beta, v_neg(beta2)), v_add(gamma, v_neg(gamma2)) )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (b2*(g-g2)*(a-a2))
    if (cur_idx <= dst_idx) and (dst_idx < nxt_idx):
        return Gilbert3D_d2xyz( dst_idx, cur_idx,
                                v_add( p, v_add( v_add( alpha, v_neg(d_alpha) ), v_add( v_add(beta2, v_neg(d_beta)), gamma2 ) ) ),
                                v_neg(beta2), v_add(gamma, v_neg(gamma2)), v_neg(v_add(alpha, v_neg(alpha2))) )
    cur_idx = nxt_idx

    return Gilbert3D_d2xyz( dst_idx, cur_idx,
                            v_add( p, v_add( v_add(alpha, v_neg(d_alpha)), v_add( gamma2, v_neg(d_gamma) ) ) ),
                            v_neg(gamma2), v_neg(v_add(alpha, v_neg(alpha2))), beta2)



def Gilbert3D_d2xyz( dst_idx, cur_idx, p, alpha, beta, gamma ):
    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)
    g = abs_sum_v(gamma)

    a0 = (a % 2)
    b0 = (b % 2)
    g0 = (g % 2)

    # base cases
    #
    if (a == 2) and (b == 2) and (g == 2):
        return Hilbert2x2x2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    if (a == 1): return Gilbert2D_d2xyz(dst_idx, cur_idx, p, beta, gamma)
    if (b == 1): return Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, gamma)
    if (g == 1): return Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, beta)

    # eccentric cases
    #
    if ((3*a) > (5*b)) and ((3*a) > (5*g)):
        return Gilbert3DS0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    if ((2*b) > (3*g)) or ((2*b) > (3*a)):
        return Gilbert3DS2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    if (2*g) > (3*b):
        return Gilbert3DS1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    # bulk recursion
    #
    if (g0 == 0):
        return Gilbert3DJ0_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    if (a0 == 0) or (b0 == 0):
        return Gilbert3DJ1_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)

    # a0 == b0 == g0 == 1
    #
    return Gilbert3DJ2_d2xyz(dst_idx, cur_idx, p, alpha, beta, gamma)


########################
#     _ ___             
#  __| |_  )_ ___  _ ___
# / _` |/ /\ \ / || |_ /
# \__,_/___/_\_\\_, /__|
#               |__/    
########################


####################################################
#   ______ ____           __  ____ ___    __    __ 
#  / ___(_) / /  ___ ____/ /_|_  // _ \__/ /___/ /_
# / (_ / / / _ \/ -_) __/ __//_ </ // /_  __/_  __/
# \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
#                                                  
####################################################


####################################################
#   ______ ____           __  ___  ___    __    __ 
#  / ___(_) / /  ___ ____/ /_|_  |/ _ \__/ /___/ /_
# / (_ / / / _ \/ -_) __/ __/ __// // /_  __/_  __/
# \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
#                                                  
####################################################

def Gilbert2D_xyz2d(cur_idx, q, p, alpha, beta):
    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)

    u = [ p[0], p[1], p[2] ]

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)

    if b == 1:
        return cur_idx + dot_v( d_alpha, v_add(q, v_neg(u)) )

    if a == 1:
        return cur_idx + dot_v( d_beta, v_add(q, v_neg(u)) )

    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)

    if (2*a) > (3*b):
        if (a2%2) and (a>2):
            alpha2 = v_add(alpha2, d_alpha)
            a2 = abs_sum_v(alpha2)

        if inBounds( q, u, alpha2, beta ):
            return Gilbert2D_xyz2d( cur_idx, q,
                                    u,
                                    alpha2,
                                    beta )
        cur_idx += a2*b
        u = v_add(u, alpha2)

        return Gilbert2D_xyz2d( cur_idx, q,
                                u,
                                v_add(alpha, v_neg(alpha2)),
                                beta )

    if (b2%2) and (b>2):
        beta2 = v_add(beta2, d_beta)
        b2 = abs_sum_v(beta2)

    if inBounds(q, u, beta2, alpha2):
        return Gilbert2D_xyz2d( cur_idx, q, u, beta2, alpha2 )
    cur_idx += b2*a2

    u = v_add(p, beta2)
    if inBounds(q, u, alpha, v_add(beta, v_neg(beta2))):
        return Gilbert2D_xyz2d( cur_idx, q,
                                u,
                                alpha,
                                v_add(beta, v_neg(beta2)) )
    cur_idx += a*(b-b2)

    u = v_add(p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_neg(d_beta)) )) 
    return Gilbert2D_xyz2d( cur_idx, q,
                            u,
                            v_neg(beta2),
                            v_add(alpha2, v_neg(alpha)) )


def Gilbert2D_d2xyz(dst_idx, cur_idx, p, alpha, beta):
    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)

    if b == 1:
        d_idx = dst_idx - cur_idx
        return v_add(p, v_mul(d_idx, d_alpha))

    if a == 1:
        d_idx = dst_idx - cur_idx
        return v_add(p, v_mul(d_idx, d_beta))

    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)

    if (2*a) > (3*b):
        if (a2%2) and (a>2):
            alpha2 = v_add(alpha2, d_alpha)
            a2 = abs_sum_v(alpha2)

        nxt_idx = cur_idx + (a2*b)
        if ((cur_idx <= dst_idx) and
            (dst_idx < nxt_idx)):
            return Gilbert2D_d2xyz( dst_idx, cur_idx, p, alpha2, beta )
        cur_idx = nxt_idx

        return Gilbert2D_d2xyz( dst_idx, cur_idx,
                                v_add(p, alpha2),
                                v_add(alpha, v_neg(alpha2)),
                                beta )

    if (b2%2) and (b>2):
        beta2 = v_add(beta2, d_beta)
        b2 = abs_sum_v(beta2)

    nxt_idx = cur_idx + (b2*a2)
    if ((cur_idx <= dst_idx) and
        (dst_idx < nxt_idx)):
        return Gilbert2D_d2xyz( dst_idx, cur_idx, p, beta2, alpha2 )
    cur_idx = nxt_idx

    nxt_idx = cur_idx + (a*(b-b2))
    if ((cur_idx <= dst_idx) and
        (dst_idx < nxt_idx)):
        return Gilbert2D_d2xyz( dst_idx, cur_idx,
                                v_add(p, beta2),
                                alpha,
                                v_add(beta, v_neg(beta2)) )
    cur_idx = nxt_idx

    return Gilbert2D_d2xyz( dst_idx, cur_idx,
                            v_add(p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_neg(d_beta)) ) ),
                            v_neg(beta2),
                            v_add(alpha2, v_neg(alpha)) )


def Gilbert2DAsync(p, alpha, beta):
    a = abs_sum_v(alpha)
    b = abs_sum_v(beta)

    d_alpha = v_delta(alpha)
    d_beta  = v_delta(beta)

    if b == 1:
        u = [ p[0], p[1], p[2] ]
        for i in range(a):
            yield u
            u = v_add(u, d_alpha)
        return

    if a == 1:
        u = [ p[0], p[1], p[2] ]
        for i in range(b):
            yield u
            u = v_add(u, d_beta)
        return

    alpha2  = v_div2(alpha)
    beta2   = v_div2(beta)

    a2 = abs_sum_v(alpha2)
    b2 = abs_sum_v(beta2)

    if (2*a) > (3*b):
        if (a2%2) and (a>2):
            alpha2 = v_add(alpha2, d_alpha)

        yield from Gilbert2DAsync( p, alpha2, beta )
        yield from Gilbert2DAsync( v_add(p, alpha2),
                                   v_add(alpha, v_neg(alpha2)),
                                   beta )
        return

    if (b2%2) and (b>2):
        beta2 = v_add(beta2, d_beta)

    yield from Gilbert2DAsync( p, beta2, alpha2 )

    yield from Gilbert2DAsync( v_add(p, beta2),
                               alpha,
                               v_add(beta, v_neg(beta2)) )

    yield from Gilbert2DAsync( v_add(p, v_add( v_add(alpha, v_neg(d_alpha)), v_add(beta2, v_neg(d_beta)) ) ),
                               v_neg(beta2),
                               v_add(alpha2, v_neg(alpha)) )
    





####################################################
#   ______ ____           __  ___  ___    __    __ 
#  / ___(_) / /  ___ ____/ /_|_  |/ _ \__/ /___/ /_
# / (_ / / / _ \/ -_) __/ __/ __// // /_  __/_  __/
# \___/_/_/_.__/\__/_/  \__/____/____/ /_/   /_/   
#                                                  
####################################################


# Note that this implementation *is* different from the original Gilbert2D
# version because of hte integer rounding.
# When coercing even/odd values, this version uses round to 0 rather
# than round to -inf, which is what the original use
#
def Gilbert2D(width, height):
    """
    Generalized Hilbert extended version ('gilbert++') space-filling curve for
    arbitrary sized 2D rectangular grids. Generates discrete 2D coordinates to fill
    a rectangle of size (width x height).
    """

    p = [0,0,0]
    alpha = [width, 0, 0]
    beta = [0, height, 0]

    yield from Gilbert2DAsync(p, alpha, beta)

def Gilbert3D(width, height, depth):
    """
    Generalized Hilbert extended version ('gilbert++') space-filling curve for
    arbitrary sized 3D rectangular grids. Generates discrete 3D coordinates to fill
    a rectangle of size (width x height x depth).
    """

    p = [0,0,0]
    alpha = [width, 0, 0]
    beta = [0, height, 0]
    gamma = [0, 0, depth]

    yield from Gilbert3DAsync(p, alpha, beta, gamma)



def spotcheck_helperfunctions():
    import random
    _irand = random.randint

    print("## sgn")
    for x in range(-3,3):
        print(x, sgn(x))
    print("")

    print("## v_delta")
    u2 = [ _irand(-10,10), _irand(-10,10) ]
    u3 = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]

    print(u2, v_delta(u2))
    print(u3, v_delta(u3))
    print("")

    print("## v_divq (3)")
    u2 = [ _irand(-10,10), _irand(-10,10) ]
    u3 = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]

    print(u2, v_divq(u2, 3))
    print(u3, v_divq(u3, 3))
    print("")

    print("## v_div2")
    u2 = [ _irand(-10,10), _irand(-10,10) ]
    u3 = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]

    print(u2, v_div2(u2))
    print(u3, v_div2(u3))
    print("")

    print("## abs_sum_v")
    u2 = [ _irand(-10,10), _irand(-10,10) ]
    u3 = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]

    print(u2, abs_sum_v(u2))
    print(u3, abs_sum_v(u3))
    print("")

    R = 9

    print("## d2e")
    for x in range(-R,R):
        print(x, d2e(x))
    print("")

    print("## d2u")
    for x in range(-R,R):
        print(x, d2u(x))
    print("")

    print("## dqe (3)")
    for x in range(-R,R):
        print(x, dqe(x, 3))
    print("")

    print("## dqu (3)")
    for x in range(-R,R):
        print(x, dqu(x, 3))
    print("")


    print("## inBounds")
    for it in range(10):
        p = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]
        q = [ _irand(-10,10), _irand(-10,10), _irand(-10,10) ]
        a = [ _irand(-100,100), 0, 0]
        b = [ 0, _irand(-100,100), 0]
        g = [ 0, 0, _irand(-100,100)]

        print(q,p,a,b,g, "inBounds:", inBounds(q,p, a,b,g))
    print("")

    print("## inBounds (2)")
    for it in range(10):
        p = [ _irand(-10,10), _irand(-10,10), 0 ]
        q = [ _irand(-10,10), _irand(-10,10), 0 ]
        a = [ _irand(-100,100), 0, 0]
        b = [ 0, _irand(-100,100), 0]

        print(q,p,a,b, "inBounds (2):", inBounds(q,p, a,b))
    print("")



if __name__ == "__main__":
    import argparse

    action = 'xy'
    depth = 1

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--action', type=str)
    parser.add_argument('width', type=int)
    parser.add_argument('height', type=int)
    parser.add_argument('-D', '--depth', type=int)
    args = parser.parse_args()

    if args.action is not None:
        action = args.action

    if args.depth is not None:
        depth = args.depth

    if action == 'xy':
        for (x,y,z) in Gilbert2D(args.width, args.height):
            print(x,y)

    elif action == 'xyz':
        for (x,y,z) in Gilbert3D(args.width, args.height, depth):
            print(x,y,z)

    elif action == 'xy2d':
        p = [0,0,0]
        alpha = [args.width,0,0]
        beta = [0,args.height,0]
        for y in range(args.height):
            for x in range(args.width):
                idx = Gilbert2D_xyz2d( 0, [x,y,0], p, alpha, beta )
                print(idx, x,y)

    elif action == 'd2xy':
        W = args.width
        H = args.height
        p = [0,0,0]
        alpha = [args.width,0,0]
        beta = [0,args.height,0]
        for idx in range(W*H):
            xyz = Gilbert2D_d2xyz( idx, 0, p, alpha, beta )
            print(xyz[0], xyz[1])

    elif action == 'd2xyz':
        W = args.width
        H = args.height
        D = depth
        p = [0,0,0]
        alpha = [args.width,0,0]
        beta = [0,args.height,0]
        gamma = [0,0,depth]
        for idx in range(W*H*D):
            xyz = Gilbert3D_d2xyz( idx, 0, p, alpha, beta, gamma )
            print(xyz[0], xyz[1], xyz[2])

    elif action == 'h2x2x2':
        alpha = [args.width,0,0]
        beta = [0,args.height,0]
        gamma = [0,0,depth]
        for (x,y,z) in Hilbert2x2x2Async([0,0,0], alpha, beta, gamma):
            print(x,y,z)

    elif action == 'test':
        spotcheck_helperfunctions()

