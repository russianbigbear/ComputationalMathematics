import matplotlib.pyplot as pl
from scipy.interpolate import CubicSpline
from scipy.misc import derivative
from numpy import arange
from numpy.random import uniform


def read_file(filename):
    """Чтение файла"""
    x, y = [], []

    # чтеник строк
    with open(filename, 'r') as f:
        lines = f.readlines()

    # деление строк на x и y
    for line in lines:
        line = line.split()
        x.append(float(line[0]))
        y.append(float(line[1]))

    return x, y


def simpson(func, l, u, auto, eps):
    """Составная формула Симпсона"""
    last = 0

    if auto:
        s = 1

    iter_count = arange(1, ((u - l) / s))

    while True:
        curr = func(l) + func(u)
        for i in iter_count:
            if i % 2 == 0:
                curr += 2 * func(l + i * s)
            else:
                curr += 4 * func(l + i * s)
        curr *= s / 3

        if abs(curr - last) < eps:
            break

        last = curr
        s /= 2
        iter_count = arange(1, ((u - l) / s))

    return curr


def monte_carlo(fx, fy, x1, x2, y1, y2, n):
    """Метод Монте-Карло"""

    def in_polygon(_x, _y, _xp, _yp):
        """входит ли в область"""
        c = 0
        for j in range(len(_xp)):
            if ((_yp[j] <= _y < _yp[j - 1]) or (_yp[j - 1] <= _y < _yp[j])) and (
                    _x > (_xp[j - 1] - _xp[j]) * (_y - _yp[j]) / (_yp[j - 1] - _yp[j]) + _xp[j]): c = 1 - c
        return c

    x = uniform(x1, x2, n)  # генерируем xi
    y = uniform(y1, y2, n)  # генерируем yi

    inside = []
    outside = []

    # проверяем принадлежность точки
    for i in range(n):
        if in_polygon(x[i], y[i], fx, fy) == 1:
            # если входит в область, то записываем в массив попаданий и переходим к следущей точке
            inside.append((x[i], y[i]))
        else:
            # если не входят в область, то записываем в масив промахов
            outside.append((x[i], y[i]))
            # S фигур * кол. попыток / общее количество
            area = abs(x1 - x2) * abs(y1 - y2) * len(inside) / n

    return area, inside, outside


def main():
    """Основная функция"""

    # вычисление функции
    def func(_x):
        return spline_x(_x) * derivative(spline_y, _x, dx=1e-10)

    # переменные
    step = 0.1  # шаг
    precision = 0.01  # точность
    points = 1000  # количество точек

    # ввод имени файла
    # красивый интерфейс file_name = input("Введите название файла: ") + ".dat"
    file_name = "input.dat"

    # чтение файла
    x, y = read_file(file_name)

    # Построеник кубического сплайна
    spline_x = CubicSpline([i for i in range(len(x))], x)
    spline_y = CubicSpline([i for i in range(len(y))], y)

    new_x = list(arange(min(x), len(x) - 1 + step, step))
    new_y = [spline_y(i) for i in new_x]

    # формирование строки + вычисление результата по формуле Симпсона
    str_simpson = "Метод Симпсона S = " + str(abs(simpson(func, min(x), len(x) - 1, True, precision)))
    t = spline_x(new_x)
    s, inside, outside = monte_carlo(t, new_y, min(t), max(t), min(new_y), max(new_y), points)

    # формирование строки  с резкльтатом по методу Монте-Карло
    str_monte_carlo = "Метод Монте-Карло S = " + str(abs(s))

    # вывод графика
    pl.plot(spline_x(new_x), new_y, color='blue', label="область")
    pl.plot(x, y, "o")

    # вывод внешних/внутрених точек
    pl.plot([i[0] for i in inside], [i[1] for i in inside], ".", color='green', label="внутрение точки")
    pl.plot([i[0] for i in outside], [i[1] for i in outside], ".", color='red' , label="внешние точки")

    # вывод строк с результатом
    pl.text(1, 14, str_simpson, bbox=dict(facecolor='pink', alpha=1))
    pl.text(1, 13, str_monte_carlo, bbox=dict(facecolor='pink', alpha=1))

    # вывод заголовка/легенды/сетки
    pl.title("10 лабораторная работа.")
    pl.legend()
    pl.grid(True)
    pl.show()


if __name__ == '__main__':
    main()
