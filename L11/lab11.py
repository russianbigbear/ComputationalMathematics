# -*- coding: utf8 -*-
# !/usr/bin/python3

STEP = 0.001  # шаг
PRECISION = 0.00001  # точность
START = -10  # начало отрезка
END = 10  # конец отрезка


def main():
    """Основная функция"""
    """Косбинированный метод, разбиение отрезка [a,b] на 3 отрезка"""

    f = eval("lambda x: x**3 - 100")  # функция
    deriv1 = eval("lambda x: 3*x*x")  # 1 производная функции
    deriv2 = eval("lambda x: 6*x")  # 2 производная функции

    a = START
    b = END

    # разбиение возможно при условии f(a)f(b) < 0
    if f(START) * f(END) < 0:
        while abs(a - b) > 2 * PRECISION:
            # Условие начальной точки для метода хорд f(x)f''(x) < 0
            if f(a) * deriv2(a) < 0:
                a -= f(a) * (a - b) / (f(a) - f(b))
            # Условие начальной точки для метода касательных f(x)f''(x) > 0
            elif f(a) * deriv2(a) > 0:
                a -= f(a) / deriv1(a)

            if f(b) * deriv2(b) < 0:
                b -= f(b) * (b - a) / (f(b) - f(a))
            elif f(b) * deriv2(b) > 0:
                b -= f(b) / deriv1(b)

    solution = (a + b) / 2

    print("Решение уравнения x^3 - 100 с заданной точностью %f = %f" % (PRECISION, solution))


if __name__ == '__main__':
    main()
