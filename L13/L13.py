# -*- coding: utf8 -*-


from numpy import *
from sympy import *

PRECISION = 0.00001  # точность
N = 4  # количество функций


def norma(x):
    """Вычисление нормы"""
    return sum([abs(x[i]) for i in range(len(x))])


def f_diff(jacob, xk):
    """Вычисление диффенциала функции"""
    return [([eval("lambda x1,x2,x3,x4: " + str(jacob[i][j]))(xk[0], xk[1], xk[2], xk[3]) for j in range(N)]) for i in range(N)]


def main():
    """Основная функция"""

    file = open('output.txt', 'w')

    # создание символьных переменных
    x1 = symbols('x1')
    x2 = symbols('x2')
    x3 = symbols('x3')
    x4 = symbols('x4')

    # предварительное обозначение функций
    f1 = "x1 - 10 * x2 + 3 * x3 * x2**2"
    f2 = "x1 + 7 * x2 + x3 - 7 * x4 + 17"
    f3 = "x1 - x2+  2 * x3 + x4 - 4"
    f4 = "x1 - x2 + x3 * x1 - x4"

    # запись в файл функций
    file.write("Функции: \n")
    file.write(f1 + "\n")
    file.write(f2 + "\n")
    file.write(f3 + "\n")
    file.write(f4 + "\n\n")

    # матрица коэфициентов (расчет Якобиана)
    jacobiyan = [[diff(f1, x1), diff(f1, x2), diff(f1, x3), diff(f1, x4)],
                 [diff(f2, x1), diff(f2, x2), diff(f2, x3), diff(f2, x4)],
                 [diff(f3, x1), diff(f3, x2), diff(f3, x3), diff(f3, x4)],
                 [diff(f4, x1), diff(f4, x2), diff(f4, x3), diff(f4, x4)]]

    # запись в файл частных производынх (коофицентов)
    file.write("Частные производные: \n")
    for item in jacobiyan:
        for i in item:
            file.write("%20s" % str(i))
        file.write("\n")
    file.write("\n")

    # система уравнений
    fe1 = eval("lambda x1,x2,x3,x4: " + f1)
    fe2 = eval("lambda x1,x2,x3,x4: " + f2)
    fe3 = eval("lambda x1,x2,x3,x4: " + f3)
    fe4 = eval("lambda x1,x2,x3,x4: " + f4)

    # создание xi
    xi = [1, 1, 1, 1]

    # коэфициенты
    fd = f_diff(jacobiyan, xi)

    # свободные члены
    f = [fe1(xi[0], xi[1], xi[2], xi[3]),
         fe2(xi[0], xi[1], xi[2], xi[3]),
         fe3(xi[0], xi[1], xi[2], xi[3]),
         fe4(xi[0], xi[1], xi[2], xi[3])]

    # решщаем СЛАУ
    xi1 = linalg.solve(fd, multiply(f, -1))

    # пока погрешность меньше точности
    while norma(xi - xi1) > PRECISION:
        xi = list(xi1)  # сохраняем предыдущее приближение
        fd = f_diff(jacobiyan, xi)  # вычисляем новые коофиценты

        # вычисляем  новые свободные члены
        f = [fe1(xi[0], xi[1], xi[2], xi[3]),
             fe2(xi[0], xi[1], xi[2], xi[3]),
             fe3(xi[0], xi[1], xi[2], xi[3]),
             fe4(xi[0], xi[1], xi[2], xi[3])]

        # решаем СЛАУ
        xi1 = linalg.solve(fd, multiply(f, -1)) + xi

    # создания списка решения
    x = list(xi1)
    # создание списка невязок
    f = [fe1(x[0], x[1], x[2], x[3]), fe2(x[0], x[1], x[2], x[3]), fe3(x[0], x[1], x[2], x[3]),
         fe4(x[0], x[1], x[2], x[3])]

    # запись решений
    file.write("Решение: \n")
    for item in x:
        file.write(str(item) + " ")
    file.write("\n\n")

    # запись невязок
    file.write("Невязки: \n")
    for item in f:
        file.write(str(item) + " ")
    file.write("\n")


if __name__ == '__main__':
    main()
