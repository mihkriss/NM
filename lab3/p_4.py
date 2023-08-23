def Diff():
    x = [0.0, 0.2, 0.4, 0.6, 0.8]
    n = len(x)
    y = [1.0, 1.4214, 1.8918, 2.4221, 3.0255]
    u = 0.4
    d = [[0 for j in range(n)] for i in range(n)]

    #Таблица разделённых разностей
    for i in range(n):
        d[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            d[i][j] = (d[i][j-1] - d[i+1][j-1]) / (x[i] - x[i+j])

    #Первая производная
    v = 0
    for j in range(2, n):
        f = 0
        for k in range(1, j):
            g = d[0][j]
            for i in range(j):
                if i != k:
                    g *= (u - x[i])
            f += g
        v += f

    #Вторая производная
    w = 0
    for j in range(3, n):
        f = 0
        for k in range((j-1)*(j-2)):
            g = d[0][j]
            for i in range(j):
                if i != 1 + (k-1)//(j-2) and i != 1 + ((k+(k-1)//(j-1)) % (j-1)):
                    g *= (u - x[i])
            f += g
        w += f
    return v, w

v, w = Diff()

print("Первая производная", round(v, 4))
print("Вторая производная", round(w, 4))
