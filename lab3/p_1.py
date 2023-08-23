import math

def lagrange_interpolation():
    v_a = 0
    v_b = 0
    u = 3*math.pi/16
    Xi_a = [math.pi/8, 2*math.pi/8, 3*math.pi/8, 4*math.pi/8]
    Xi_b = [math.pi/8, math.pi/3, 3*math.pi/8, math.pi/2]
    n = len(Xi_a)
    
    y_a = [1/math.tan(x)+x for x in Xi_a]
    y_b = [1/math.tan(x)+x for x in Xi_b]
    w_a = [0]*n
    w_b = [0]*n

    for j in range(n):
        w_a[j] = y_a[j]
        f = y_a[j]
        for i in range(n):
            if i != j:
                w_a[j] = w_a[j]*(u-Xi_a[i])/(Xi_a[j]-Xi_a[i]) # интерполяционные многочлены Лагранжа
                f = f*(u-Xi_a[i])/(Xi_a[j]-Xi_a[i])
        v_a += f
        y_p = u*math.log(u)
        L_a = abs(v_a-y_p) # значение погрешности интерполяции в точке X*

    for j in range(n):
        w_b[j] = y_b[j]
        f = y_b[j]
        for i in range(n):
            if i != j:
                w_b[j] = w_b[j]*(u-Xi_b[i])/(Xi_b[j]-Xi_b[i]) # интерполяционные многочлены Лагранжа
                f = f*(u-Xi_b[i])/(Xi_b[j]-Xi_b[i])
        v_b += f
        y_p = u*math.log(u)
        L_b = abs(v_b-y_p) # значение погрешности интерполяции в точке X*
    
    return w_a, L_a, w_b, L_b


w_a, L_a, w_b, L_b = lagrange_interpolation()

print("Вычисления коэффициентов интерполяционного многочлена Лагранжа: ")
print('---------------------------------------')
print('Для а): ')
print('---------------------------------------')
print('Интерполяционные многочлены Лагранжа: ')
print([round(num, 4) for num in w_a])
print('Значение погрешности интерполяции в точке X*: ')
print(round(L_a, 4))
print('---------------------------------------')

print('Для б): ')
print('---------------------------------------')
print('Интерполяционные многочлены Лагранжа: ')
print([round(num, 4) for num in w_b])
print('Значение погрешности интерполяции в точке X*: ')
print(round(L_b, 4))
print('---------------------------------------')


def Newton():
    v_a = 0
    v_b = 0
    u = 3*math.pi/16
    Xi_a = [math.pi/8, 2*math.pi/8, 3*math.pi/8, 4*math.pi/8]
    Xi_b = [math.pi/8, math.pi/3, 3*math.pi/8, math.pi/2]
    n = len(Xi_a)
    
    y_a = [1/math.tan(x)+x for x in Xi_a]
    y_b = [1/math.tan(x)+x for x in Xi_b]
    w_a = [0]*n
    w_b = [0]*n
    d_a = [[0]*n for i in range(n)]
    d_b = [[0]*n for i in range(n)]
     
    for i in range(n):
        d_a[i][0] = y_a[i]
    for j in range(1, n):
        for i in range(n-j):
            d_a[i][j] = (d_a[i][j-1] - d_a[i+1][j-1]) / (Xi_a[i] - Xi_a[i+j])
    for j in range(n):
        w_a[j] = d_a[0][j]
        f = d_a[0][j]
        for i in range(j):
            f *= (u - Xi_a[i])
        v_a += f
        y_p = u*math.log(u)
        L_a = abs(v_a - y_p) # значение погрешности интерполяции в точке u
   
    for i in range(n):
        d_b[i][0] = y_b[i]
    for j in range(1, n):
        for i in range(n-j):
            d_b[i][j] = (d_b[i][j-1] - d_b[i+1][j-1]) / (Xi_b[i] - Xi_b[i+j])
    for j in range(n):
        w_b[j] = d_b[0][j]
        f = d_b[0][j]
        for i in range(j):
            f *= (u - Xi_b[i])
        v_b += f
        y_p = u*math.log(u)
        L_b = abs(v_b - y_p) # значение погрешности интерполяции в точке u
    return w_a, L_a, w_b, L_b

w_a, L_a, w_b, L_b = Newton()

print("Вычисления коэффициентов интерполяционного многочлена Ньютона: ")
print('---------------------------------------')
print('Для а): ')
print('---------------------------------------')
print('Интерполяционные многочлены Ньютона: ')
print([round(num, 4) for num in w_a])
print('Значение погрешности интерполяции в точке X*: ')
print(round(L_a, 4))
print('---------------------------------------')

print('Для б): ')
print('---------------------------------------')
print('Интерполяционные многочлены Ньютона: ')
print([round(num, 4) for num in w_b])
print('Значение погрешности интерполяции в точке X*: ')
print(round(L_b, 4))
print('---------------------------------------')
