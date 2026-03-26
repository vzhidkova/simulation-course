import math
import numpy as np


x0 = 0.4
m = 1.0
a = 5.6
N = 10 ** 5

my_data = []
my_data.append(x0)
library_data = np.random.uniform(0.0, m, N)
M1 = x0
M2 = 0

for i in range(N):
    x0 = (x0 * a) % m
    my_data.append(x0)
    M1 += x0
    M2 += library_data[i]

M1 = M1/N
M2 = M2/N
D1 = D2 = 0
for i in range(N):
    D1 += (my_data[i] - M1) ** 2
    D2 += (library_data[i] - M2) ** 2

D1 = D1/N
D2 = D2/N


print("Теоретические результаты:", m/2, m**2/12)
print("Мультипликативный конгруэнтный генератор:", M1, D1)
print("Встроенный генератор:", M2, D2)

