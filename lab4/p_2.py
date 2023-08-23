import math
import numpy as np
import matplotlib.pyplot as plt

PI = math.pi

def func(x, y, y_der):
    return 2 * (1 + (math.tan(x) ** 2)) * y


def g(x, y, k):
    return k


def p(x):
    return 0


def q(x):
    return - 2 * (1 + (math.tan(x) ** 2))


def pureFunction(x):
    return - math.tan(x)


def f(x):
    return 0


def first(x, y, x0):
    i = 0
    while i < len(x) - 1 and x[i + 1] < x0:
        i += 1
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i])

def rungeKutt(f, a, b, h, y0, y_der):
    n = int((b - a) / h)
    x = [i for i in np.arange(a, b + h, h)]
    y = [y0]
    k = [y_der]
    for i in range(n):
        K1 = h * g(x[i], y[i], k[i])
        L1 = h * f(x[i], y[i], k[i])
        K2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * K1, k[i] + 0.5 * L1)
        L2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K1, k[i] + 0.5 * L1)
        K3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * K2, k[i] + 0.5 * L2)
        L3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K2, k[i] + 0.5 * L2)
        K4 = h * g(x[i] + h, y[i] + K3, k[i] + L3)
        L4 = h * f(x[i] + h, y[i] + K3, k[i] + L3)
        y.append(y[i] + (K1 + 2 * K2 + 2 * K3 + K4) / 6)
        k.append(k[i] + (L1 + 2 * L2 + 2 * L3 + L4) / 6)
    return x, y, k

def stop(y, y1, eps):
    if abs(y[-1] - y1) > eps:
        return True
    else:
        return False


def newN(n_last, n, ans_last, ans, b, y1):
    x, y = ans_last[0], ans_last[1]
    phi_last = y[-1] - y1
    x, y = ans[0], ans[1]
    phi = y[-1] - y1
    return n - (n - n_last) / (phi - phi_last) * phi


def shootingMethod(a, b, y0, y1, h, eps):
    n_last = 1
    n = 0.8
    y_der = n_last
    ans_last = rungeKutt(func, a, b, h, n_last, y_der)[:2]
    y_der = n
    ans = rungeKutt(func, a, b, h, n, y_der)[:2]

    while stop(ans[1], y1, eps):
        n, n_last = newN(n_last, n, ans_last, ans, b, y1), n
        ans_last = ans
        y_der = n
        ans = rungeKutt(func, a, b, h, y0, y_der)[:2]

    return ans


def tma(a, b, c, d, shape):
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    x = [0] * (shape + 1)
    for i in range(1, shape):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    for i in reversed(range(shape)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x[:-1]


def finiteDifferenceMethod(a, b, alpha, beta, delta, gamma, y0, y1, h):
    n = int((b - a) / h)
    x = [i for i in np.arange(a, b + h, h)]
    A = [0] + [1 - p(x[i]) * h / 2 for i in range(0, n - 1)] + [-gamma]
    B = [alpha * h - beta] + [q(x[i]) * h ** 2 - 2 for i in range(0, n - 1)] + [delta * h + gamma]

    C = [beta] + [1 + p(x[i]) * h / 2 for i in range(0, n - 1)] + [0]
    D = [y0 * h] + [f(x[i]) * h ** 2 for i in range(0, n - 1)] + [y1 * h]

    y = tma(A, B, C, D, len(A))
    return x, y



def rungeRomberg(ans, exact):
    k = ans[0]['h'] / ans[1]['h']
    Y1 = [yi for xi, yi in zip(ans[0]['Shooting']['x'], ans[0]['Shooting']['y']) if xi in ans[1]['Shooting']['x']]
    Y2 = [yi for xi, yi in zip(ans[1]['Shooting']['x'], ans[1]['Shooting']['y']) if xi in ans[0]['Shooting']['x']]
    shoot_err = [y1 + (y2 - y1) / (k ** 2 - 1) for y1, y2 in zip(Y1, Y2)]
    X_ex = [xi for xi in ans[0]['Shooting']['x'] if xi in ans[1]['Shooting']['x']]
    Y_ex = [pureFunction(i) for i in X_ex]
    for i in range(len(shoot_err)):
        shoot_err[i] = abs(shoot_err[i] - Y_ex[i])

    Y1 = [yi for xi, yi in zip(ans[0]['FD']['x'], ans[0]['FD']['y']) if xi in ans[1]['FD']['x']]
    Y2 = [yi for xi, yi in zip(ans[1]['FD']['x'], ans[1]['FD']['y']) if xi in ans[0]['FD']['x']]
    fd_err = [y1 + (y2 - y1) / (k ** 2 - 1) for y1, y2 in zip(Y1, Y2)]
    X_ex = [xi for xi in ans[0]['FD']['x'] if xi in ans[1]['FD']['x']]
    Y_ex = [pureFunction(i) for i in X_ex]
    for i in range(len(fd_err)):
        fd_err[i] = abs(fd_err[i] - Y_ex[i])

    return {'Shooting': shoot_err, 'FD': fd_err}


def sse(f, y):
    return round(sum([(f_i - y_i) ** 2 for f_i, y_i in zip(f, y)]), 5)


if __name__ == '__main__':

    a = 0
    b = PI / 6
    alpha = 1
    delta = 1
    gamma = 0
    beta = 0
    y0 = 0.0
    y1 = pureFunction(b)
    step = PI / 30
    eps = 1e-5

    print(f'Interval: [{a}, {b}]')
    print(f'y0 = {y0}, y1 = {y1}')
    print()

    res = []
    res2 = []
    ans = []
    steps = [step, step / 2]
    i = 0

    for h in steps:
        print(f'Current step: {h}')
        print('Shooting method')
        res.append(shootingMethod(a, b, y0, y1, h, eps))
        for x, y in zip(res[i][0], res[i][1]):
            print(f'x: {round(x, 5)}, y: {round(y, 5)}')
        print()

        print('Finite difference method')
        res2.append(finiteDifferenceMethod(a, b, alpha, beta, delta, gamma, y0, y1, h))
        for x, y in zip(res2[i][0], res2[i][1]):
            print(f'x: {round(x, 5)}, y: {round(y, 5)}')
        print()

        ans.append({
            "h": h,
            "Shooting": {'x': res[i][0],  'y': res[i][1]},
            "FD":       {'x': res2[i][0], 'y': res2[i][1]}
        })

        i += 1

    exact = []
    for h in steps:
        x_ex = [i for i in np.arange(a, b + h, h)]
        y_ex = [pureFunction(i) for i in x_ex]
        exact.append((x_ex, y_ex))

    err = rungeRomberg(ans, exact)
    print("All errors")
    print('Shooting method runge error: {}'.format(err['Shooting']))
    print('Finite difference method runge error: {}'.format(err['FD']))
  