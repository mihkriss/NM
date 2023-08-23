import math

# Функция правой части ОДУ
def f(x, y, z):
    return -z/x

# Точное решение ОДУ
def y_exact(x):
    return 1 + math.log(abs(x))

# Метод Эйлера
def euler_method(x0, y0, z0, h, n):
    x = x0
    y = y0
    z = z0
    for i in range(n):
        z += h*f(x, y, z)
        y += h*z
        x += h
    return y

# Метод Рунге-Кутты 4-го порядка
def runge_kutta_method(x0, y0, z0, h, n):
    x = x0
    y = y0
    z = z0
    for i in range(n):
        k1 = h*z
        l1 = h*f(x, y, z)
        k2 = h*(z + 0.5*l1)
        l2 = h*f(x + 0.5*h, y + 0.5*k1, z + 0.5*l1)
        k3 = h*(z + 0.5*l2)
        l3 = h*f(x + 0.5*h, y + 0.5*k2, z + 0.5*l2)
        k4 = h*(z + l3)
        l4 = h*f(x + h, y + k3, z + l3)
        z += (k1 + 2*k2 + 2*k3 + k4)/6
        y += (l1 + 2*l2 + 2*l3 + l4)/6
        x += h
    return y

# Метод Адамса 4-го порядка
def adams_method(x0, y0, z0, h, n):
    x = x0
    y = y0
    z = z0
    y_prev = y
    z_prev = z
    for i in range(3):
        k1 = h*z
        l1 = h*f(x, y, z)
        k2 = h*(z + 0.5*l1)
        l2 = h*f(x + 0.5*h, y + 0.5*k1, z + 0.5*l1)
        k3 = h*(z + 0.5*l2)
        l3 = h*f(x + 0.5*h, y + 0.5*k2, z + 0.5*l2)
        k4 = h*(z + l3)
        l4 = h*f(x + h, y + k3, z + l3)
        z_next = z + (k1 + 2*k2 + 2*k3 + k4)/6
        y_next = y + (l1 + 2*l2 + 2*l3 + l4)/6
        x += h
        y_prev = y
        z_prev = z
        y = y_next
        z = z_next
    return y

#Вычисление погрешности метода Рунге-Ромберга
def runge_romberg(h, y1, y2, p):
    return (y2 - y1)/(2**p - 1)

#Задание начальных условий
x0 = 1
y0 = 1
z0 = 1
h = 0.1
a = 1
b = 2
#Число шагов сетки
n = int((b - a)/0.1)

#Метод Эйлера
y_euler = euler_method(x0, y0, z0, h, n)
print("Метод Эйлера: y(2) = %.4f" % y_euler)
err = runge_romberg(h, y_exact(2), y_euler, 1)
print("Погрешность метода Рунге-Ромберга: %.4f" % err)
print('----------------------------------------')

#Метод Рунге-Кутты 4-го порядка
y_rk4 = runge_kutta_method(x0, y0, z0, h, n)
print("Метод Рунге-Кутты 4-го порядка: y(2) = %.4f" % y_rk4)
err = runge_romberg(h, y_exact(2), y_rk4, 4)
print("Погрешность метода Рунге-Ромберга: %.4f" % err)
print('----------------------------------------')

#Метод Адамса 4-го порядка
y_adams = adams_method(x0, y0, z0, h, n)
print("Метод Адамса 4-го порядка: y(2) = %.4f" % y_adams)
err = runge_romberg(h, y_exact(2), y_adams, 4)
print("Погрешность метода Рунге-Ромберга: %.4f" % err)
print('----------------------------------------')

y_exact_2 = y_exact(2)
print("Точное решение: y(2) = %.4f" % y_exact_2)
print('----------------------------------------')

#Сравнение с точным решением
print("Погрешности относительно точного решения:")
print("Метод Эйлера: %.4f" % abs(y_exact_2 - y_euler))
print("Метод Рунге-Кутты 4-го порядка: %.4f" % abs(y_exact_2 - y_rk4))
print("Метод Адамса 4-го порядка: %.4f" % abs(y_exact_2 - y_adams))

