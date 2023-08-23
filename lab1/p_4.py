import math

def Jacobi(A, eps):
    n = len(A)  # размерность матрицы

    # Инициализация матриц вращения
    V = [[int(i == j) for j in range(n)] for i in range(n)]
    a = [[A[i][j] for j in range(n)] for i in range(n)]

    # Вычисление максимального недиагональный элемента матрицы:
    def max_element():
        max_el = 0
        l, m = 1, 2
        for i in range(n):
            for j in range(i+1, n):
                if abs(a[i][j]) > max_el:
                    max_el = abs(a[i][j])
                    l, m = i, j
        return max_el, l, m

    # Основной цикл
    k = 0
    while True:
        max_el, l, m = max_element()
        if max_el < eps:
            break

        # Вычисление угла вращения
        phi = 0.5 * math.atan(2 * a[l][m] / (a[l][l] - a[m][m]))

        # Вычисление матрицы вращения
        U = [[int(i == j) for j in range(n)] for i in range(n)]
        U[l][l] = math.cos(phi)
        U[l][m] = -math.sin(phi)
        U[m][l] = math.sin(phi)
        U[m][m] = math.cos(phi)

        # Обновление матриц
        b = [[0 for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                b[i][j] = 0
                for p in range(n):
                    b[i][j] += U[p][i] * a[p][j]
        for i in range(n):
            for j in range(n):
                a[i][j] = 0
                for p in range(n):
                    a[i][j] += b[i][p] * U[p][j]

        Vt = [[U[j][i] for j in range(n)] for i in range(n)]
        b = [[0 for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                for p in range(n):
                    b[i][j] += V[i][p] * Vt[p][j]
        V = b
        
        k += 1
    
    # Получение собственных значений и собственных векторов
    w = [a[i][i] for i in range(n)]
    v = [[V[j][i] for j in range(n)] for i in range(n)]

    return w, v, k




A = [[-3, 3, -1], [-1, 8, 1], [3, 1, 5]]
eps = 0.2
w, v, k = Jacobi(A, eps)

for i in range(len(v)):
    for j in range(len(v[i])):
        v[i][j] = round(v[i][j], 3) 
print("Собственные значения:")
print([round(num, 3) for num in w])
print("Собственные векторы:")
print(v)
print("Количество итераций:" )
print(k)
