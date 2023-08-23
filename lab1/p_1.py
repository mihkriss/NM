import numpy as np
import matplotlib.pyplot as plt

def u(x, y, a, b, N):
    result = 0
    for n in range(1, N+1):
        B_n = 1 / np.sinh(n*np.pi*b/a)
        result += (B_n * np.sin(n*np.pi*x/a) * np.sinh(n*np.pi*y/a)) / np.sinh(n*np.pi*b/a)
    return result

a = 1  # Длина прямоугольника
b = 1  # Ширина прямоугольника
N = 20  # Количество гармоник для суммирования

x = np.linspace(0, a, 100)
y = np.linspace(0, b, 100)
X, Y = np.meshgrid(x, y)
Z = u(X, Y, a, b, N)

plt.contourf(X, Y, Z, cmap='viridis')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('График функции u(x, y)')
plt.show()
