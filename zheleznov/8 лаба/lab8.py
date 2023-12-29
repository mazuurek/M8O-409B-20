import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from functions import *
import math
import logging

a = 1
T = 3
n = 15
k = 10
m = 10
hx = math.pi/(4*n)
hy = math.log(2)/m
tau = T/k

t_array = [0] * (k+1)
x_array = [0] * (n+1)
y_array = [0] * (m+1)
solution_array_t = [0] * (k+1)

for i in range(k+1):
    t_array[i] = tau * i
for i in range(n+1):
    x_array[i] = hx * i
for i in range(m+1):
    y_array[i] = hy * i
def T_y_back(t,y):
    return math.cosh(y) * math.exp(-3 * a * t)
def T_y_front(t,y):
    return 0
def X_t_left(x,t):
    return math.cos(2 * x) * math.exp(-3 * a * t)
def X_t_right(x,t):
    return (3/4) * math.cos(2 * x) * math.exp(-3 * a * t)
def X_y_bottom(x,y):
    return math.cos(2 * x) * math.cosh(y)
def solution(x,y,t):
    return math.cos(2 * x) * math.cosh(y) * math.exp(-3 * a * t)

# сетка
array = [0] * (n+1)
for i in range(n+1):
    array[i] = [0]*(m+1)
    for j in range(m+1):
        array[i][j] = [0] * (k+1)


for i in range(m+1):
    for j in range(n+1):
        array[j][i][0] = X_y_bottom(j * hx, i * hy)

for i in range(m+1):
    for j in range(k+1):
        array[0][i][j] = T_y_back(tau * j, hy * i)

for i in range(n+1):
    for j in range(k+1):
        array[i][0][j] = X_t_left(i * hx, j * tau)

for i in range(n+1):
    for j in range(k+1):
        array[i][-1][j] = X_t_right(i * hx, j * tau)

#Переменные направления

for iter in range(1,k+1):
    A = [0] * (n - 1)
    for i in range(n - 1):
        A[i] = [0] * (n - 1)

    C = [0] * (m - 1)
    for i in range(m - 1):
        C[i] = [0] * (m - 1)

    for i in range(n - 1):
        A[i][i] = 1 + ((a * tau) / (hx**2))
    count = 0
    for i in range(1, n - 1):
        A[i][count] = - (a * tau) / (2 * hx**2)
        count += 1
    count = 1
    for i in range(n - 2):
        A[i][count] = - (a * tau) / (2 * hx**2)
        count += 1

    for i in range(m - 1):
        C[i][i] = 1 + (a * tau) / (hy**2)
    count = 0
    for i in range(1, m - 1):
        C[i][count] = - (a * tau) / (2 * hy**2)
        count += 1
    count = 1
    for i in range(m - 2):
        C[i][count] = - (a * tau) / (2 * hy**2)
        count += 1

    semi_array = [0] * (n + 1)
    for i in range(n + 1):
        semi_array[i] = [0] * (m + 1)

    for i in range(n+1):
        semi_array[i][0] = X_t_left(hx * i,tau * ((iter - 1) + 0.5)) #iter с 1 начинается,проверь!
    for i in range(m+1):
        semi_array[0][i] = T_y_back(tau * ((iter - 1) + 0.5),i * hy)
    for i in range(n+1):
        semi_array[i][-1] = X_t_right(hx * i,tau * ((iter - 1) + 0.5))

    B = [0] * (n - 1)
    D = [0] * (m - 1)

    for j in range(1,m-1):
        for i in range(1,n-2):
            B[i] = array[i+1][j][iter-1] + ((a * tau) / (2 * hy**2)) * (array[i+1][j-1][iter-1] - (2 * array[i+1][j][iter-1]) + array[i+1][j+1][iter-1])
        B[0] = ((a * tau) / (2 * hx**2)) * semi_array[0][j] + array[1][j][iter-1] + ((a * tau) / (2 * hy**2)) * (array[1][j-1][iter-1] - (2 * array[1][j][iter-1]) + array[1][j+1][iter-1])
        B[-1] = ((a * tau) / (2 * hx**2)) * semi_array[-1][j] + array[n-1][j][iter-1] + ((a * tau) / (2 * hy**2)) * (array[n-1][j-1][iter-1] - (2 * array[n-1][j][iter-1]) + array[n-1][j+1][iter-1])
        x = progonka(A, B, n-1)
        for i in range(len(x)):
            semi_array[i+1][j] = x[i]

    for j in range(1,n-1):
        for i in range(1,m-2):
            D[i] = semi_array[j][i+1] + ((a * tau) / (2 * hx**2)) * (semi_array[j+1][i+1] - 2 * semi_array[j][i+1] + semi_array[j-1][i+1])
        D[0] = ((a * tau) / (2 * hy**2)) * array[j][0][iter] + semi_array[j][1] + ((a * tau) / (2 * hx**2)) * (semi_array[j+1][1] - 2 * semi_array[j][1] + semi_array[j-1][1])
        D[-1] = ((a * tau) / (2 * hy**2)) * array[j][-1][iter] + semi_array[j][m-1] + ((a * tau) / (2 * hx**2)) * (semi_array[j+1][m-1] - 2 * semi_array[j][m-1] + semi_array[j-1][m-1])
        y = progonka(C, D, m-1)
        for i in range(len(y)):
            array[j][i+1][iter] = y[i]


ax = plt.axes()
ax.set_xlabel("t")
ax.set_ylabel("u(x,y,t)")
for i in range(k+1):
    solution_array_t[i] = solution(x_array[5],y_array[6],t_array[i])
print(solution_array_t)
plt.plot(t_array, solution_array_t, label="Аналитическое решение")
plt.plot(t_array, array[5][6], label=f"Вычисленное решение (k={k})")
plt.legend(loc='best')
plt.title("Метод переменных направлений")
plt.show()


#дробные шаги

for iter in range(1,k+1):
    A = [0] * (n - 1)
    for i in range(n - 1):
        A[i] = [0] * (n - 1)

    C = [0] * (m - 1)
    for i in range(m - 1):
        C[i] = [0] * (m - 1)

    for i in range(n - 1):
        A[i][i] = 1 + ((2 * a * tau) / (hx**2))
    count = 0
    for i in range(1, n - 1):
        A[i][count] = - (a * tau) / (hx**2)
        count += 1
    count = 1
    for i in range(n - 2):
        A[i][count] = - (a * tau) / (hx**2)
        count += 1

    for i in range(m - 1):
        C[i][i] = 1 + (2 * a * tau) / (hy**2)
    count = 0
    for i in range(1, m - 1):
        C[i][count] = - (a * tau) / (hy**2)
        count += 1
    count = 1
    for i in range(m - 2):
        C[i][count] = - (a * tau) / (hy**2)
        count += 1

    semi_array = [0] * (n + 1)
    for i in range(n + 1):
        semi_array[i] = [0] * (m + 1)

    for i in range(n+1):
        semi_array[i][0] = X_t_left(hx * i,tau * ((iter - 1) + 0.5)) #iter с 1 начинается,проверь!
    for i in range(m+1):
        semi_array[0][i] = T_y_back(tau * ((iter - 1) + 0.5),i * hy)
    for i in range(n+1):
        semi_array[i][-1] = X_t_right(hx * i,tau * ((iter - 1) + 0.5))

    B = [0] * (n - 1)
    D = [0] * (m - 1)

    for j in range(1,m-1):
        for i in range(1,n-2):
            B[i] = array[i+1][j][iter-1]
        B[0] = ((a * tau) / (hx**2)) * semi_array[0][j] + array[1][j][iter-1]
        B[-1] = ((a * tau) / (hx**2)) * semi_array[-1][j] + array[n-1][j][iter-1]
        x = progonka(A, B, n-1)
        for i in range(len(x)):
            semi_array[i+1][j] = x[i]

    for j in range(1,n-1):
        for i in range(1,m-2):
            D[i] = semi_array[j][i+1]
        D[0] = ((a * tau) / (hy**2)) * array[j][0][iter] + semi_array[j][1]
        D[-1] = ((a * tau) / (hy**2)) * array[j][-1][iter] + semi_array[j][m-1]
        y = progonka(C, D, m-1)
        for i in range(len(y)):
            array[j][i+1][iter] = y[i]

ax = plt.axes()
ax.set_xlabel("t")
ax.set_ylabel("u(x,y,t)")
for i in range(k+1):
    solution_array_t[i] = solution(x_array[5],y_array[6],t_array[i])
print(solution_array_t)
plt.plot(t_array, solution_array_t, label="Аналитическое решение")
plt.plot(t_array, array[5][6], label=f"Вычисленное решение (k={k})")
plt.legend(loc='best')
plt.title("Метод дробных шагов")
plt.show()