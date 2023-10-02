import numpy as np
from numpy import *


def all_finite(z1, z2, z3, w1, w2, w3):
    a = np.linalg.det(
        [[z1*w1, w1, 1],
         [z2*w2, w2, 1],
         [z3*w3, w3, 1]])
    b = np.linalg.det(
        [[z1*w1, z1, w1],
         [z2*w2, z2, w2],
         [z3*w3, z3, w3]])
    c = np.linalg.det(
        [[z1, w1, 1],
         [z2, w2, 1],
         [z3, w3, 1]])
    d = np.linalg.det(
        [[z1*w1, z1, 1],
         [z2*w2, z2, 1],
         [z3*w3, z3, 1]])
    print("\nall_finite\n")
    print("all_finite a,b,c,d = ", a, b, c, d, "\n")
    return (a, b, c, d)


def z1_infinite(z1, z2, z3, w1, w2, w3):
    assert (np.isinf(z1))
    a = w1*(w2 - w3)
    b = w1*(z2*w3 - z3*w2) + w2*w3*(z3 - z2)
    c = w2 - w3
    d = w1*(z2 - z3) - z2*w2 + z3*w3
    print("\nz1_infinite\n")
    print("z1_infinite a,b,c,d = ", a, b, c, d, "\n")
    return (a, b, c, d)


def w1_infinite(z1, z2, z3, w1, w2, w3):
    assert (np.isinf(w1))
    a = z1*(w2 - w3) - z2*w2 + z3*w3
    b = z1*(z2*w3 - z3*w2) + z2*z3*(w2 - w3)
    c = z3 - z2
    d = z1*(z2 - z3)
    print("\nw1_infinite\n")
    print("w1_infinite a,b,c,d = ", a, b, c, d, "\n")
    return (a, b, c, d)


def z1w1_infinite(z1, z2, z3, w1, w2, w3):
    assert (np.isinf(z1) and np.isinf(w1))
    a = w2 - w3
    b = z2*w3 - z3*w2
    c = 0
    d = z2 - z3
    print("\nz1w1_infinite\n")
    print("z1w1_infinite a,b,c,d = ", a, b, c, d, "\n")
    return (a, b, c, d)


def z1w2_infinite(z1, z2, z3, w1, w2, w3):
    assert (np.isinf(z1) and np.isinf(w2))
    a = w1
    b = -z2*w3 + z3*(w3 - w1)
    c = 1
    d = -z2
    print("\nz1w2_infinite\n")
    print("z1w2_infinite a,b,c,d = ", a, b, c, d, "\n")
    return (a, b, c, d)


def mobius_coeff(z1, z2, z3, w1, w2, w3):

    infz = np.isinf(z1) or np.isinf(z2) or np.isinf(z3)
    infw = np.isinf(w1) or np.isinf(w2) or np.isinf(w3)

    if infz:
        if np.isinf(z2):
            z1, z2 = z2, z1
            w1, w2 = w2, w1
        if np.isinf(z3):
            z1, z3 = z3, z1
            w1, w3 = w3, w1
        if infw:
            if np.isinf(w1):
                return z1w1_infinite(z1, z2, z3, w1, w2, w3)
            if np.isinf(w3):
                z2, z3 = z3, z2
                w2, w3 = w3, w2
            return z1w2_infinite(z1, z2, z3, w1, w2, w3)
        else:
            return z1_infinite(z1, z2, z3, w1, w2, w3)

    if infw:  # and all z finite
        if np.isinf(w2):
            z1, z2 = z2, z1
            w1, w2 = w2, w1
        if np.isinf(w3):
            z1, z3 = z3, z1
            w1, w3 = w3, w1
        return w1_infinite(z1, z2, z3, w1, w2, w3)

    return all_finite(z1, z2, z3, w1, w2, w3)


def mobius(x, a, b, c, d):
    if np.isinf(x):
        if c == 0:
            return np.inf
        return a/c
    if c*x + d == 0:
        return np.inf
    else:
        return (a*x + b)/(c*x + d)


def test_mobius(z1, z2, z3, w1, w2, w3):
    tolerance = 1e-6
    a, b, c, d = mobius_coeff(z1, z2, z3, w1, w2, w3)
    print("z1,z2,z3,w1,w2,w3 :", z1, z2, z3, w1, w2, w3, "\n")
    print("mobius_coeff: ", a, b, c, d, "\n")
    for (x, y) in [(z1, w1), (z2, w2), (z3, w3)]:
        m = mobius(x, a, b, c, d)
        assert (np.isinf(m) and np.isinf(y) or abs(m - y) <= tolerance)


# test_mobius(np.inf, 1, 0, 2, 1/2 + 3j/2, 2j-1)
test_mobius(1, 1j, -1j, 1j+1, 0, 2)
