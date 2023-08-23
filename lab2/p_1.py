import math
import numpy as np
import matplotlib.pyplot as plt

def Iter():
    x = x0
    k = 0
    while True:
        x_new = phi(x)
        k += 1
        if abs(x_new - x) * q/(1 - q) < e:
            return x_new, k
        x = x_new

def phi(x):
    return math.sqrt(2-2**x)

e = 0.0001
q = 0.98
a = 0.6
b = 0.75
x0 = (a + b)/2

root, k = Iter()
print('------------------------------------------')
print("Метод простых итераций: ")
print("Корень: ", round(root, 4))
print("Количество итераций: ", k)
print('------------------------------------------')

def Newton():
    x = x0
    k = 0
    if (df(a))*(ddf(a)) > 0:
        x = a
    else:
        if (df(b))*(ddf(b)) > 0:
            x = b
    while True:
        x_new = x - f(x)/ df(x)
        k += 1
        if abs(x_new - x) < e:
            return x_new, k
        x = x_new

def f(x):
    return 2**x + x**2 - 2
def df(x):
    return (2**x)*math.log(2) + 2*x
def ddf(x):
    return 2*math.log(2) + 1 

e = 0.0001
a = 0.6
b = 0.75
x0 = (a + b)/2

root, iterations = Newton()
print("Метод Ньютона: ")
print("Корень: ", round(root, 4))
print("Количество итераций: ", iterations)
print('------------------------------------------')


#График
def func(x):
    return 2**x + x**2 - 2

plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams.update({'font.size': 18})

x = np.linspace(-2, 3, 100)

plt.plot(x, func(x), label=r'2^x + x^2 - 2')

plt.axhline(y=0, color='black', lw=2)
plt.axvline(x=0, color='black', lw=2)

# добавляем вертикальную линию в точке корня
root = 0.6535
plt.axvline(x=root, color='red', linestyle='--')

# добавляем горизонтальную линию, пересекающую вертикальную линию
plt.axhline(y=0.25, color='green', linestyle='--')

# добавляем точку на пересечении линий
plt.plot(root, 0.25, 'bo', markersize=10)

# подписываем оси
plt.xlabel('x')
plt.ylabel('y')

plt.legend()
plt.show()