#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:11:59 2024

@author: ben
"""

from numpy import array, matrix, zeros, diag, diagflat, dot, matmul
from numpy.linalg import norm, inv

def jacobi(
        A: array, 
        b: array, 
        x: array, 
        e: float) -> array:
    """
    Solves the equation Ax=b using the Jacobi iterative method.

    Parameters:
    -----------
    A : array
        Coefficient matrix.
    b : array
        Right-hand side of the system of equations.
    x : array
        Initial vector.
    e : float
        Error condition.

    Returns:
    --------
    x : array
        Solution vector.
    """                                                                                                                                                                  
    D = diag(A)
    R = A - diagflat(D)
    error = norm(dot(A, x) - b)

    while error > e:
        x = dot(inv(diagflat(D)), (b - dot(R, x)))
        error = norm(dot(A, x) - b)

    return x

def gaussSeidel(
        A: array, 
        b: array, 
        x: array, 
        e: float) -> array:
    """
    Solves the equation Ax=b using the Gauss-Seidel method.

    Parameters:
    -----------
    A : array
        Coefficient matrix.
    b : array
        Right-hand side of the system of equations.
    x : array
        Initial vector.
    e : float
        Error condition.

    Returns:
    --------
    x : array
        Solution vector.
    """
    error = norm(dot(A, x) - b)
    n = b.shape[0]
    xnew = zeros((n, 1))

    while error > e:
        x = xnew.copy()
        for i in range(0, n):
            xnew[i] = (b[i] - 
                       sum([A[i, j] * x[j] for j in range(i + 1, n)]) - 
                       sum([A[i, j] * xnew[j] for j in range(i)])) / A[i, i]
        error = norm(dot(A, x) - b)

    return x


if __name__ == "__main__":
    A = array([
        [4, 1, 0],
        [1, 4, 1],
        [0, 1, 4]
    ])
    b = array([[1], [1], [1]])
    x0 = array([[1], [2], [3]])
    e = 1e-6

    # Solve using Jacobi method
    x = jacobi(A, b, x0, e)
    print("Jacobi solution:\n", x)

    # Solve using Gauss-Seidel method
    x = gaussSeidel(A, b, x0, e)
    print("Gauss-Seidel solution:\n", x)