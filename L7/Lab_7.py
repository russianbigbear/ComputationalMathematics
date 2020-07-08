# !/usr/bin/python3
import math
from collections import defaultdict
import matplotlib.pyplot as plt
from numpy import arange
import bisect


def e(x):
    """Вычисление y для экспоненты"""
    return math.exp(x)


def sinus(x):
    """Вычисление y для синуса"""
    return  math.sin(x)


def sinus4(x):
    """Вычисление y  для синуса*4"""
    return  math.sin(4*x)


def xxx(x):
    """Вычисление y  для x^3"""
    return  x * x * x

def atans(x):
    """Вычисление y"""
    return  math.atan(50*x)


# Класс "точка"
class Dot:
    def __init__(self, x, y): self.x, self.y = [x, y]


# Структура, описывающая сплайн на каждом сегменте сетки
class Tuple:
    a, b, c, d, x = [0., 0., 0., 0., 0.]

splines = defaultdict(lambda: Tuple())

def build_spline(dots):
    """
    Построение слайна
    В dots:
        x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
        y - значения функции в узлах сетки
    """

    # инициализация структуры сплайнов
    for i in range(len(dots)):
        splines[i].x, splines[i].a = dots[i].x, dots[i].y
    alpha, beta = [defaultdict(lambda: 0.), defaultdict(lambda: 0.)]

    # решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
    # вычисление прогоночных коэффициентов - прямой ход метода прогонки
    for i in range(1, len(dots) - 1):
        hi = dots[i].x - dots[i-1].x
        hi1 = dots[i+1].x - dots[i].x
        A = hi
        C = 2. * (hi + hi1)
        B = hi1
        F = 6. * ((dots[i + 1].y - dots[i].y) / hi1 - (dots[i].y - dots[i - 1].y) / hi)
        z = (A * alpha[i-1] + C)
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z

    # нахождение решения - обратный ход метода прогонки
    for i in reversed(range(1, len(dots) - 1)):
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i]

    # по известным коэффициентам c[i] находим значения b[i] и d[i]
    for i in reversed(range(1, len(dots))):
        hi = dots[i].x - dots[i - 1].x  # h[i] = x[i] - x[i-1]
        splines[i].d = (splines[i].c - splines[i - 1].c) / hi  # d[i] = (c[i] - c[i-1])/h[i]
        splines[i].b = (hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 +
                        (dots[i].y - dots[i - 1].y) / hi)  # b[i] = h[i]*c[i]/2 - (h[i]^2)*d[i]/6 +  (f[i] - f[i-1])/h[i]
    return splines


def lagrange(x_list, y_list, x):
    """Вычисление полинома Лагранжа в заданном отрезка, с заданными точками."""
    if len(x_list) != len(y_list):
        raise ValueError("Списки координат не равны!!!")

    y = 0
    for j in range(len(y_list)):
        p1, p2 = 1, 1
        #считаем произведение(П) (x-xi)/(xj-xi)
        for i in range(len(x_list)):
            if i != j:
                p1 *= x - x_list[i]
                p2 *= x_list[j] - x_list[i]
        y += y_list[j] * p1 / p2  # возвращаем сумму yi * на П
    return y


def interpolate(x):
    """Вычисление значения интерполированной функции в произвольной точке"""
    distribution = sorted([t[1].x for t in splines.items()])
    indx = bisect.bisect_left(distribution, x)
    if indx == len(distribution):
        return 0
    dx = x - splines[indx].x
    return (splines[indx].a + splines[indx].b * dx +
            splines[indx].c * dx ** 2 / 2 +
            splines[indx].d * dx ** 3 / 6)

def main():
    """Основная функция"""
    file_name = "input.dat"

    # выбор функции для вычисления y по задданому x
    func = {"e": e, "sin": sinus, "sin4": sinus4, "xxx": xxx, "atan": atans}
    func_name = input("Введите название функции (e, sin, sin4, xxx, atan): ")

    # чтение из файла
    with open(file_name) as file:
        n = int(file.readline())
        x, y = [], []

        temp = [float(num) for num in  file.readline().split()]
        for i in range(n):
            x.append(temp[i])
            y.append(func[func_name](x[i]))

    # упорядочивание элементов по значениям и создание списка
    temp = sorted(zip(x, y), key=lambda c: c[0])
    x, y = [i[0] for i in temp], [i[1] for i in temp]

    #построение сплайна
    splines = build_spline([Dot(x[i], y[i]) for i in range(n)])

    # создание списка x от x[min] до x[max]
    new_x = list(arange(min(x), max(x), 0.1))
    new_x.append(new_x[-1] + 0.1)

    new_y = [func[func_name](x) for x in new_x]

    new_y_lag = [lagrange(x, y, i) for i in new_x]
    new_y_inter = [interpolate(x) for x in new_x]

    # вывод коофицентов a,b,c,d,x для каждой точки
    print("%35s" % "Коофиценты:")
    print("%10s %10s %10s %10s %10s " % ("a", "b", "c", "d", "x" ))
    [print("%10f %10f %10f %10f %10f " % (splines[i].a, splines[i].b, splines[i].c, splines[i].d, splines[i].x) ) for i in splines]

    #Вывод графика
    plt.plot(x, y, "o", label="точки")
    plt.plot(new_x, new_y, "b", label = "функция")
    plt.plot(new_x, new_y_lag, "r", label = "Лагранж")
    plt.plot(new_x, new_y_inter, "g", label = "интерполяция")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()

