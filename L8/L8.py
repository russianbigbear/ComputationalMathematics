from sympy import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

def max_in_list(lst, s, e):
    """Максимальное в списке"""
    assert lst
    m = abs(lst[s])
    for i in range(s, len(lst)-e):
        if abs(lst[i]) > m:
            m = abs(lst[i])
    return m


def derivative1(x, y, step):
    """Вычисление первой производной"""

    # u'(x0) = (-3u(x0) + 4u(x1) - u(x2)) / 2h
    res = [(-3 * y[0] + 4 * y[1] - y[2]) / (2 * step)]

    # u'(xi) = (u(x[i+1]) - u(x[i-1])) / 2h
    for i in range(1, len(y) - 1):
        res.append((y[i + 1] - y[i - 1]) / (2 * step))

    return res + [(y[-1] - y[-2]) / (x[-1] - x[-2])] # результат вычислений


def derivative2(y, h):
    """Вычисление второй производной"""
    res = [float(0)]

    # u''(xi) = (u(x[i-1] - 2u(xi) + u(x[i+1]) / h^2
    for i in range(1, len(y) - 1):
        res.append((y[i - 1] - 2 * y[i] + y[i + 1]) / (h ** 2))

    return res + [float(0)] # результат вычислений


def derivative3(y, h):
    """Вычисление третьей производной"""
    res = [float(0), float(0)]

    h = h**3 * 2 # для быстроты вычисление h = 2h^3

    # u'''(xi) = (2u(x[i-1] - u(x[i-2]) - 2u(x[i+1]) + u(x[i+2])) / 2h^3
    for i in range(2, len(y) - 2):
        res.append((2*y[i - 1] - y[i - 2] - 2*y[i + 1] + y[i + 2]) / h)

    return res + [float(0), float(0)] # результат вычислений


def main():
    """Основная функция"""
    h = 0.5 # шаг в таблице значений
    x_min = 0  # минимальное значение аргумента
    x_max = 5  # максимальное значение аргумента

    x = symbols('x')
    func = "sin(x)"
    print("Функция: " + func)
    differential1 = diff(func, x)
    print("Первая производная: " + str(differential1))
    differential2 = diff(differential1, x)
    print("Вторая производная: " + str(differential2))
    differential3 = diff(differential2, x)
    print("Третья производная: " + str(differential3))

    f = eval("lambda x:" + str(func))
    d1 = eval("lambda x:" + str(differential1))
    d2 = eval("lambda x:" + str(differential2))
    d3 = eval("lambda x:" + str(differential3))

    arguments = []
    y = []
    dif1 = []
    dif2 = []
    dif3 = []
    temp = x_min

    # заполнение данными
    while temp < x_max:
        arguments.append(temp)  # x
        y.append(f(temp))  # y = sin(x)
        dif1.append(d1(temp))   # u'
        dif2.append(d2(temp))  # u''
        dif3.append(d3(temp))  # u'''
        temp += h

    # вычисление производных по таблицам
    derivat1 = derivative1(arguments, y, h)
    derivat2 = derivative2(y, h)
    derivat3 = derivative3(y, h)

    # вывод результатов в файл(1)
    with open("output1.dat", 'w') as file:
        file.write('%20s' % 'Точка')
        file.write('%20s' % 'Функция')
        file.write('%20s' % 'Произ. аналит. 1')
        file.write('%20s' % 'Произ. аналит. 2')
        file.write('%20s' % 'Произ. аналит. 3')
        file.write('%20s' % 'Произ. числ. 1')
        file.write('%20s' % 'Произ. числ. 2')
        file.write('%20s' % 'Произ. числ. 3')
        file.write('\n')

        for i in range(len(arguments)):
            file.write("%20.14f" % arguments[i])
            file.write("%20.14f" % y[i])
            file.write("%20.14f" % dif1[i])
            file.write("%20.14f" % dif2[i])
            file.write("%20.14f" % dif3[i])
            file.write("%20.14f" % derivat1[i])
            file.write("%20.14f" % derivat2[i])
            file.write("%20.14f\n" % derivat3[i])

    # отрисовка графиков
    plt.plot(arguments, y, "r", label = "Функция")
    plt.plot(arguments, dif1, "b", label = "Произ. анал. 1")
    plt.plot(arguments, derivat1, "g", label = "Произ. числ. 1")
    plt.grid(True)
    plt.xlabel("Ось x")
    plt.ylabel("Ось y")
    plt.title("Первый график")
    plt.legend()
    plt.show()

    plt.plot(arguments, y, "r", label = "Функция")
    plt.plot(arguments, dif2, "b", label = "Произ. анал. 2")
    plt.plot(arguments, derivat2, "g", label = "Произ. числ. 2")
    plt.grid(True)
    plt.xlabel("Ось x")
    plt.ylabel("Ось y")
    plt.title("Второй график")
    plt.legend()
    plt.show()


    plt.plot(arguments, y, "r", label = "Функция")
    plt.plot(arguments, dif3, "b", label = "Произ. анал. 3")
    plt.plot(arguments, derivat3, "g", label = "Произ. числ. 3")
    plt.grid(True)
    plt.xlabel("Ось x")
    plt.ylabel("Ось y")
    plt.title("Третий график")
    plt.legend()
    plt.show()

    # вывод результатов в файл (2)
    with open("output2.dat", 'w') as file:
        file.write('%20s' % 'Точка')
        file.write('%20s' % 'Функция')
        file.write('%20s' % 'Разность по 1')
        file.write('%20s' % 'Разность по 2')
        file.write('%20s' % 'Разность по 3')
        file.write('\n')

        norm1 = [abs(dif1[i] - derivat1[i]) for i in range(len(arguments))]
        max_in_list(norm1, 1, len(norm1)-1)
        norm2 = [abs(dif2[i] - derivat2[i]) for i in range(len(arguments))]
        max_in_list(norm2, 2, len(norm1)-2)
        norm3 = [abs(dif3[i] - derivat3[i]) for i in range(len(arguments))]
        max_in_list(norm3, 3, len(norm1)-3)

        for i in range(len(arguments)):
            file.write("%20.14f" % arguments[i])
            file.write("%20.14f" % y[i])
            file.write("%20.14f" % norm1[i])
            file.write("%20.14f" % norm2[i])
            file.write("%20.14f" % norm3[i])
            file.write('\n')

    # Внести в таблицу значений возмущения, не превышающие delta по модулю
    delta = 0.001

    for i in range(len(y)):
        y[i] = y[i] + (np.random.randint(100)-50) * delta

    interfer_d1 = derivative1(arguments, y, h)
    interfer_d2 = derivative2(y, h)
    interfer_d3 = derivative3(y, h)

    # вывод результатов в файл (3)
    with open("output3.dat", 'w') as file:
        file.write('%20s' % 'Точка')
        file.write('%25s' % 'Функция с возмущением')
        file.write('%25s' % 'Разность по 1')
        file.write('%20s' % 'Разность по 2')
        file.write('%20s' % 'Разность по 3')
        file.write('\n')

        for i in range(len(arguments)):
            file.write("%20.14f" % arguments[i])
            file.write("%25.14f" % y[i])
            file.write("%25.14f" % abs(derivat1[i] - interfer_d1[i]))
            file.write("%20.14f" % abs(derivat2[i] - interfer_d2[i]))
            file.write("%20.14f" % abs(derivat3[i] - interfer_d3[i]))
            file.write('\n')

if __name__ == '__main__':
    main()