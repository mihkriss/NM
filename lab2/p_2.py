import math
import numpy as np
import matplotlib.pyplot as plt

def Simple_itor():
    a = 0.5
    b = 0.75
    c = 0.5
    d = 0.75
    e = 0.0001
    k = 0
    q = 0.005
    r = 2*e

    x = (a+b)/2
    y = (c+d)/2

    def phi_1(x, y):
        return math.cos(y) + 3
    def phi_2(x,y):
        return math.sin(x) + 3

    while r > e:
        u = x
        v = y

        # Для метода простой итерации
        x = phi_1(u, v)
        y = phi_2(u, v)

        r = max(abs(x - u), abs(y - v))
        r = r * q / (1 - q)
        k += 1

    return x, y, k

x, y, k = Simple_itor()
print('------------------------------------------')
print("Метод простых итераций: ")
print('------------------------------------------')
print("Корень по x1: ", round(x, 4))
print("Корень по x2: ", round(y, 4))
print("Количество итераций: ", k)
print('------------------------------------------')

def Newton():
    a = 0.5
    b = 0.75
    c = 0.5
    d = 0.75
    e = 0.0001
    k = 0
    r = 2*e

    x = a+b/2
    y = c+d/2

    def f(x,y):
        return x - math.cos(y) - 3
    
    def g(x, y):
        return y - math.sin(x) - 3

    def dfdx(x, y):
        return 1
    
    def dfdy(x, y):
        return math.sin(y)
    
    def dgdx(x, y):
        return -math.cos(x)
    
    def dgdy(x, y):
        return 1

    while r > e:
        u = x
        v = y
        x = u - (f(u,v)*dgdy(u,v)-g(u,v)*dfdy(u,v))/(dfdx(u,v)*dgdy(u,v)-dfdy(u,v)*dgdx(u,v))
        y = v - (dfdx(u,v)*g(u,v)-f(u,v)*dgdx(u,v))/(dfdx(u,v)*dgdy(u,v)-dfdy(u,v)*dgdx(u,v))
        r = max(abs(x - u), abs(y - v))
        k += 1

    return x, y, k

x, y, k = Newton()
print("Метод Ньютона: ")
print('------------------------------------------')
print("Корень по x1: ", round(x, 4))
print("Корень по x2: ", round(y, 4))
print("Количество итераций: ", k)
print('------------------------------------------')

#График
x = np.linspace(-2, 2, 1000)
y1 = np.cos(x) + 3
y2 = np.sin(x) + 3
plt.xlabel('x1')
plt.ylabel('x2')
plt.plot(x, y1, label='f(x,y)')
plt.plot(x, y2, label='g(x,y)')
plt.legend()
plt.grid()

root_x = 3.729
plt.axhline(y=root_x, color='green', linestyle='--')
root_y = 0.785
plt.axvline(x=root_y, color='red', linestyle='--')

plt.axhline(y=0, color='black', lw=2)
plt.axvline(x=0, color='black', lw=2)

# Находим координаты точки пересечения
idx = np.argwhere(np.diff(np.sign(y1 - y2))).flatten()
x_intersect = x[idx]
y_intersect = y1[idx]

# Отмечаем точку на графике
plt.plot(x_intersect, y_intersect, 'ro')
plt.show()










    



