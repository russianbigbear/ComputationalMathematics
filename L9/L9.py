from sympy import *
from math import *
import matplotlib.pyplot as pl
import numpy as np

MAX_SPLIT = 200


def method_of_rectangles(func, min_lim, max_lim, delta, init, auto=True):
    """Метод прямоугольников"""

    # вычисление интеграла по формуле средних прямоугольников
    def integrate_func(_function, _min_lim, _max_lim, _n):
        integral_value = 0.0  # начальное значение интеграла
        step_integrate = (_max_lim - _min_lim) / _n  # шаг

        # разбили интервал на части с шагом step_integrate
        for _x in np.arange(_min_lim, _max_lim, step_integrate):
            # формула средних прямоугольников интеграл[a,b](f(x)dx) = h * sum(f(xi - h/2)
            integral_value += step_integrate * _function(_x + step_integrate / 2)
        return integral_value

    integral = 0.0  # начальное значение интеграла
    total = init  # число отрезков в разбиении
    step = (max_lim - min_lim) / total  # шаг

    # интегрирвание с автоматическим выбором шага
    if auto is True:
        # разбиваем интервал на xi с определенным шагом
        for x in np.arange(min_lim, max_lim, step):
            d, n = 1, 1
            # цикл с опреденной точностью и количеством разбиений
            while fabs(d) > delta and n < MAX_SPLIT:
                total += n # подсчет разбиений
                # Правило Рунге
                d = (integrate_func(func, x, x + step, n * 2) - integrate_func(func, x, x + step, n)) / 3
                n *= 2
            integral += fabs(integrate_func(func, x, x + step, n)) + d
    # интегрирование с постояным шагом
    else:
        # разбиваем интервал на xi с определенным шагом
        for x in np.arange(min_lim, max_lim, step):
            # Правило Рунге об ошибке в оценке определенного интеграла (1/3)|I2n - In|
            d = (integrate_func(func, x, x + step, 2) - integrate_func(func, x, x + step, 1)) / 3
            # Вычисление интеграла по методу
            integral += fabs(integrate_func(func, x, x + step, 1)) + d

    # вывод
    print("Метод прямоугольников: %d шага(ов), интеграл = %.5f" % (total, integral))


def trapezium_method(func, min_lim, max_lim, delta, init, auto=True):
    """Метод трапеций"""

    # вычисление интеграла по формуле
    def integrate_func(_function, _min_lim, _max_lim, _n):
        integral_value = 0.0
        step_integrate = (_max_lim - _min_lim) / _n

        # разбили интервал на части с шагом step_integrate
        for _x in np.arange(_min_lim, _max_lim, step_integrate):
            # формула трапеций интеграл[a,b](f(x)dx) = h( (f(x0) + f(xn) / 2)
            integral_value += step_integrate * (_function(_x) + _function(_x + step_integrate)) / 2
        return integral_value

    # переменные
    integral = 0.0
    total = init
    step = (max_lim - min_lim) / total

    # интегрирвание с автоматическим выбором шага
    if auto is True:
        for x in np.arange(min_lim, max_lim, step):
            d, n = 1, 1
            while fabs(d) > delta and n < MAX_SPLIT:
                total += n
                d = (integrate_func(func, x, x + step, n * 2) - integrate_func(func, x, x + step, n)) / 3
                n *= 2
            integral += fabs(integrate_func(func, x, x + step, n)) + d
    # интегрирование с постояным шагом
    else:
        for x in np.arange(min_lim, max_lim, step):
            # Правило Рунге об ошибке в оценке определенного интеграла (1/3)|I2n - In|
            d = (integrate_func(func, x, x + step, 2) - integrate_func(func, x, x + step, 1)) / 3
            # Вычисление интеграла по методу
            integral += fabs(integrate_func(func, x, x + step, 1)) + d

    print("Метод трапеций: %d шага(ов), интеграл = %.5f" % (total, integral))


def simpson_method(func, min_lim, max_lim, delta, init, auto=True):
   """Метод Симпсона"""

   # вычисление интеграла по формуле
   def integrate_func(_function, _min_lim, _max_lim, _n):
       integral_value = 0.0
       step_integrate = (_max_lim - _min_lim) / _n

       # разбили интервал на части с шагом step_integrate
       for _x in np.arange(_min_lim + step_integrate / 2, _max_lim + step_integrate / 2, step_integrate):
           h =  step_integrate / 6
           fi_m1 = _function(_x - step_integrate / 2)  # f(xi-1)
           fi_p1 = _function(_x + step_integrate / 2)  # f(xi+1)
           fi = 4*_function(_x)  # 4 * f(xi)

           # формула трапеций интеграл[a,b](f(x)dx) = h/6(f(xi-1) + 4f(xi) + f(xi+1))
           integral_value += h * (fi_m1 + fi + fi_p1)
       return integral_value

   # переменные
   integral = 0.0
   total = init
   step = (max_lim - min_lim) / total

   # интегрирвание с автоматическим выбором шага
   if auto is True:
       for x in np.arange(min_lim, max_lim , step):
           d, n = 1, 1
           while fabs(d) > delta and n < MAX_SPLIT:
               total += n
               d = (integrate_func(func, min_lim, max_lim, n * 2) - integrate_func(func, min_lim, max_lim, n)) / 15
               n *= 2
           integral += fabs(integrate_func(func, x, x + step, n)) + d
   # интегрирование с постояным шагом
   else:
       for x in np.arange(min_lim, max_lim, step):
           # Правило Рунге об ошибке в оценке определенного интеграла (1/15)|I2n - In|
           d = (integrate_func(func, x, x + step, 2) - integrate_func(func, x, x + step, 1)) / 15
           # Вычисление интеграла по методу
           integral += fabs(integrate_func(func, x, x + step, 1)) + d

   print("Метод Симпсона: %d шага(ов), интеграл = %.5f" % (total, integral))


def main():
    """Основная функция"""
    # минимальное значение аргумента
    x_min = float(input("Введите нижнюю границу: "))
    # максимальное значение аргумента
    x_max = float(input("Введите верхнюю границу: "))
    # число отрезков в разбиении
    split = 300
    # точность
    precision = 0.0001
    func = str(input("Введите функцию: "))
    f = eval("lambda x:" + func)

    print("Заданная функция: " + func)
    print("Пределы интегрирования: от %.3f до %.3f\n" % (x_min, x_max))

    """Вычисление интеграла 3 методами"""
    method_of_rectangles(f, x_min, x_max, precision, split, True)  # метод прямоугольников
    trapezium_method(f, x_min, x_max, precision, split, True)  # метод трапеций
    simpson_method(f, x_min, x_max, precision, split, True)  # метод Симпсона

    # вывод графика
    pl.plot([x for x in np.arange(x_min, x_max, precision)], [f(x) for x in np.arange(x_min, x_max, precision)], "r", label = func)
    pl.grid(true)
    pl.legend()
    pl.show()


if __name__ == '__main__':
    main()
