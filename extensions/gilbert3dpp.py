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


######
#
# helper
# functions
#
######

def sgn(t): return 1 if (t>0) else (-1 if (t<0) else 0)

# vector element wise absolute value
#
def v_delta(v): return list(map( lambda t: sgn(t), v ))
def v_divq(v,q): return list(map( lambda t: int(t/q), v ))
def v_div2(v): return v_divq(v,2)

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
    if (T2%2) == 0: return s*T2
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

######
#
# helper
# functions
#
######


def Gilbert2dAsync(x, y, ax, ay, bx, by):
    pass

if __name__ == "__main__":
    pass

