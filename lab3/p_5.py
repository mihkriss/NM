from scipy.integrate import quad
import numpy as np

# Определяем функцию, которую мы будем интегрировать
def f(x):
    return x/(x**4+81)

# Задаем границы интегрирования
X0 = 0
Xk = 2

# Вычисляем точное значение интеграла с помощью функции quad
I_exact, _ = quad(f, X0, Xk)

# Метод прямоугольников
def rectangle_rule(f, X0, Xk, n):
    h = (Xk - X0) / n
    x = np.linspace(X0 + h/2, Xk - h/2, n)
    return h * np.sum(f(x))

# Метод трапеций
def trapezoidal_rule(f, X0, Xk, n):
    h = (Xk - X0) / n
    x = np.linspace(X0, Xk, n+1)
    return h/2 * (np.sum(f(x)) - f(X0) - f(Xk))

# Метод Симпсона 
def simpsons_rule(f, X0, Xk, n):
    h = (Xk - X0) / n
    x = np.linspace(X0, Xk, n+1)
    y = f(x)
    return h/3 * (np.sum(y[0:-1:2]) + 4*np.sum(y[1::2]) + y[-1])

# Задаем значения шага
h1 = 0.5
h2 = 0.25

# Вычисляем интегралы методами прямоугольников, трапеций и Симпсона
I_rectangle_h1 = rectangle_rule(f, X0, Xk, int((Xk-X0)/h1))
I_rectangle_h2 = rectangle_rule(f, X0, Xk, int((Xk-X0)/h2))
I_trapezoidal_h1 = trapezoidal_rule(f, X0, Xk, int((Xk-X0)/h1))
I_trapezoidal_h2 = trapezoidal_rule(f, X0, Xk, int((Xk-X0)/h2))
I_simpsons_h1 = simpsons_rule(f, X0, Xk, int((Xk-X0)/h1))
I_simpsons_h2 = simpsons_rule(f, X0, Xk, int((Xk-X0)/h2))

#Вычисляем оценку погрешности с помощью Метода Рунге-Ромберга
p = 2 # порядок метода
I_RR_rectangle = I_rectangle_h1 + ((I_rectangle_h1 - I_rectangle_h2) / (2**p - 1)) - I_exact
I_RR_trapezoidal = I_trapezoidal_h1 + ((I_trapezoidal_h1 - I_trapezoidal_h2) / (2**p - 1)) - I_exact
I_RR_simpsons = I_simpsons_h1 + ((I_simpsons_h1 - I_simpsons_h2) / (2**p - 1)) - I_exact

#Вычисляем абсолютную погрешность интерполяции для каждого метода
err_rectangle = np.abs(I_exact - I_RR_rectangle)
err_trapezoidal = np.abs(I_exact - I_RR_trapezoidal)
err_simpsons = np.abs(I_exact - I_RR_simpsons)

print('Точное значение интеграла: {}'.format(round(I_exact, 7)))
print('-------------------------------------------------------------')
print('Метод прямоугольников с h1 = {}: {}'.format(round(h1, 7), round(I_rectangle_h1, 7)))
print('Метод прямоугольников с h2 = {}: {}'.format(round(h2, 7), round(I_rectangle_h2, 7)))
print('Абсолютную погрешность интерполяции для прямоугольника: {}'.format(round(err_rectangle, 7)))
print('Оценка погрешности метода прямоугольников: {}'.format(round(I_RR_rectangle, 7)))
print('-------------------------------------------------------------')
print('Метод трапеций с h1 = {}: {}'.format(round(h1, 7), round(I_trapezoidal_h1, 7)))
print('Метод трапеций с h2 = {}: {}'.format(round(h2, 7), round(I_trapezoidal_h2, 7)))
print('Абсолютную погрешность интерполяции для трапеции: {}'.format(round(err_trapezoidal, 7)))
print('Оценка погрешности метода трапеций: {}'.format(round(I_RR_trapezoidal, 7)))
print('-------------------------------------------------------------')
print('Метод Симпсона с h1 = {}: {}'.format(round(h1, 7), round(I_simpsons_h1, 7)))
print('Метод Симпсона с h2 = {}: {}'.format(round(h2, 7), round(I_simpsons_h2, 7)))
print('Абсолютную погрешность интерполяции для Симпсона: {}'.format(round(err_simpsons, 7)))
print('Оценка погрешности метода Симпсона: {}'.format(round(I_RR_simpsons, 7)))
