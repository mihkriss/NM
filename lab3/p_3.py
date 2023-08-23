import numpy as np
import matplotlib.pyplot as plt

# Таблично заданная функция
x = np.array([1.0, 1.9, 2.8, 3.7, 4.6, 5.5])
y = np.array([3.4142, 2.9818, 3.3095, 3.8184, 4.3599, 4.8318])


# Приближающий многочлен 1-ой степени
A = np.vstack([x, np.ones(len(x))]).T
a1, a0 = np.linalg.lstsq(A, y, rcond=None)[0]
P1 = a1*x + a0
F1 = np.sum((y - P1)**2)

# Приближающий многочлен 2-ой степени
A = np.vstack([x**2, x, np.ones(len(x))]).T
a2, a1, a0 = np.linalg.lstsq(A, y, rcond=None)[0]
P2 = a2*pow(x, 2) + a1*x + a0
F2 = np.sum((y - P2)**2)

print('Приближающий многочлен 1-ой степени: ')
print([round(num, 4) for num in P1])
print(round(F1, 4))

print('Приближающий многочлен 2-ой степени: ')
print([round(num, 4) for num in P2])
print(round(F2, 4))

# Создание массива x_new и вычисление соответствующих значений P1 и P2
x_new = np.linspace(0, 6, 1000)
P2_new = a2*pow(x_new, 2) + a1*x_new + a0

# Построение графиков
fig, ax = plt.subplots()
ax.plot(x, y, 'o', label='Таблично заданная функция')
ax.plot(x, P1, label=f'Приближающий многочлен 1-ой степени\nF1 = {F1:.4f}')
ax.plot(x_new, P2_new, label=f'Приближающий многочлен 2-ой степени\nF2 = {F2:.4f}')
ax.legend()
plt.show()
