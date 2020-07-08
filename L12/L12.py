# -*- coding: utf8 -*-
# !/usr/bin/python3


import math


PRECISION = 0.00001  # погрешность
START = 1.  # начало отрезка
END = 5.  # конец отрезка


def f(x):
    """Вычисление функции"""
    return (x * x) - math.log(x) - (2 * math.cos(x)) - (math.tan(x) * math.sin(x * x * x))


def _c(a, b):
    """Вычисление середины отрезка"""
    return (a + b) / 2.0  # середина отрезка


def _A(a, b, c):
    """Расчитываем коофицент А"""
    return (((f(b) - f(c)) / (b - c)) - (f(c) - f(a)) / (c - a)) / (b - a)


def _B(a, b, c):
    """Расчитываем коофицент В"""
    return (f(c) - f(a)) / (c - a) + ((((f(b) - f(c)) / (b - c)) - ((f(c) - f(a)) / (c - a))) / (b - a)) * (a - c)


def _C(a):
    """Расчитываем коофицент С"""
    return f(a)


def find_x():
    """нахождение x"""
    a = START  # левая граница интервала
    b = END  # правая граница интервала
    c = _c(a, b)  # подсчет середины интервала

    k = 0  # шаг итерации
    xi = [0]
    x = 0.0  # корень

    eps = abs(a - b)  # начальная точнсть

    if f(a) * f(b) < 0:
        # итерации проводятся до тех пор, пока |xn-1  -  xn| < eps
        while eps >= PRECISION:
            k += 1
            maximum = math.sqrt(_B(a, b, c) * _B(a, b, c) - 4 * _A(a, b, c) * _C(a))
            # нахождения 2ух точек пересечение Ох интерполяционным многочленом P2(x)
            x1 = a - (2 * _C(a) / (_B(a, b, c) + maximum))
            x2 = a - (2 * _C(a) / (_B(a, b, c) - maximum))

            # выбираем значение находящиеся в интервале
            if a < x1 < b:
                x = x1
            else:
                x = x2

            # задание новых границ для точного нахождения приближенного
            if f(a) * f(x) < 0:
                b = x
            elif f(x) * f(b) < 0:
                a = x

            xi.append(x)
            eps = abs(xi[k] - xi[k - 1])  # новая точность

    return x, k, xi


def main():
    """Основная функция"""

    result, step, xi = find_x()
    print("На %d шаге итерации приближенное значение корня нелинейного уравнения будет определяться значением %f" % (step, result))


if __name__ == '__main__':
    main()
