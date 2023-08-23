def Cubic_spline(): 
    x = [1.0, 1.9, 2.8, 3.7, 4.6]
    y = [2.8069, 1.8279, 1.6091, 1.5713, 1.5663]
    n = len(x)
    u = 2.66666667
    a = [0]*(n)
    b = [0]*(n)
    c = [0]*(n)
    d = [0]*(n)
    p = [0]*(n-1)
    q = [0]*(n-1)
    for i in range(n):
        a[i] = y[i]
    c[0] = 0

    # Метод прогонки для коэффициентов C
    p[0] = -(x[2]-x[1])/(2*x[2]-2*x[0])
    q[0] = 3*((y[2]-y[1])/(x[2]-x[1])-(y[1]-y[0])/(x[1]-x[0]))/(2*x[2]-2*x[0])

    for i in range(1, n-2):
        p[i] = -(x[i+2]-x[i+1])/((2*x[i+2]-2*x[i])+(x[i+1]-x[i])*p[i-1])
        q[i] = (3*((y[i+2]-y[i+1])/(x[i+2]-x[i+1])-(y[i+1]-y[i])/(x[i+1]-x[i]))-(x[i+1]-x[i])*q[i-1])/((2*x[i+2]-2*x[i])+(x[i+1]-x[i])*p[i-1])
    c[n-1] = q[n-2]

    for i in range(n-3, 0, -1):
        c[i]=p[i]*c[i+1]+q[i]
    
    for i in range(n-1):
        b[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(x[i+1]-x[i])*(c[i+1]+2*c[i])/3
        d[i] = (c[i+1]-c[i])/(x[i+1]-x[i])/3

    b[n-1] = (y[n-1]-y[n-2])/(x[n-1]-x[n-2])-(x[n-1]-x[n-2])*2*c[n-2]/3
    d[n-1] = -c[n-2]/(x[n-1]-x[n-2])/3

    j = 0
    for i in range(n-1):
        if u >= x[i] and u <= x[i+1]:
            j = i
    u -= x[j]
    v = a[j]+b[j]*u+c[j]*u*u+d[j]*u*u*u
    return a, b, c, d, v

a, b, c, d, v = Cubic_spline()

print('_______________________________________')
print('Кубический сплайн: ')
print('---------------------------------------')
print("a_i: ")
print([round(num, 4) for num in a])
print('---------------------------------------')
print("b_i: ")
print([round(num, 4) for num in b])
print('---------------------------------------')
print("c_i: ")
print([round(num, 4) for num in c])
print('---------------------------------------')
print("d_i: ")
print([round(num, 4) for num in d])
print('---------------------------------------')
print("Значение функции в точке X*: ")
print(round(v, 4))
print('_______________________________________')
