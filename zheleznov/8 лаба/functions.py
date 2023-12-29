import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import math

def column(mtr,i):
    return [row[i] for row in mtr]
def print_net(array,k,n):
    for i in range(k + 1):
        for j in range(n + 1):
            print(round(array[i][j], 3), " ", end="")
        print("\n")
    print()
def difference(solution,x_array,time_array,array,n,k):
    sum = 0.0
    for i in range(k+1):
        for j in range(n+1):
            sum += (solution(x_array[j],time_array[i])-array[i][j])**2
    print(f"Среднеквадратичное отклоненние:{sum/((n+1)*(k+1))}")

def progonka(A, B,n):
    P= [0] * n
    Q= [0] * n
    P[0]=-A[0][1]/A[0][0]
    Q[0]=B[0]/A[0][0]

    for i in range(1,n-1):
        P[i]=-A[i][i+1]/(A[i][i]+A[i][i-1]*P[i-1])
        Q[i]=(B[i]-A[i][i-1]*Q[i-1])/(A[i][i]+A[i][i-1]*P[i-1])

    P[n-1]=0
    Q[n-1]=(B[n-1]-A[n-1][n-2]*Q[n-2])/(A[n-1][n-1]+A[n-1][n-2]*P[n-2])

    x = [0] * n
    x[n-1]=Q[n-1]

    for i in range(n-2,-1,-1):
        x[i]=P[i]*x[i+1]+Q[i]
    return x
def init(array,U_0,U_left,U_right,n,k,h,tau):
    for i in range(n+1): # было n
        array[0][i] = U_0(i*h)

    for i in range(k):
        array[i+1][0] = U_left((i+1)*tau)

    for i in range(k):
        array[i+1][n] = U_right((i+1)*tau) # вот тут вот
def print_solution(array,time_array,x_array,solution,k,scheme_name):
    ax = plt.axes()
    ax.set_xlabel("t")
    ax.set_ylabel("u(x,t)")
    solution_array_t = [0] * (k + 1)
    plt.plot(time_array, column(array, 6), label=f"Вычисленное решение (k={k})")
    for i in range(k + 1):
        solution_array_t[i] = solution(x_array[6], time_array[i])
    print(solution_array_t)
    plt.plot(time_array, solution_array_t, label="Аналитическое решение")
    plt.legend(loc='best')
    plt.title(f"{scheme_name}")
    plt.show()
# def print_solution_x(array,time_array,x_array,solution,n,scheme_name):
#     ax = plt.axes()
#     ax.set_xlabel("x")
#     ax.set_ylabel("u(x,t)")
#     solution_array_t = [0] * (n + 1)
#     plt.plot(x_array, array[16], label="Вычисленное решение")
#     for i in range(n + 1):
#         solution_array_t[i] = solution(x_array[i], time_array[16])
#
#     plt.plot(x_array, solution_array_t, label="Аналитическое решение")
#     plt.legend(loc='best')
#     plt.title(f"{scheme_name}")
#     plt.show()