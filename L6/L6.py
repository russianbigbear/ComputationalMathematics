import numpy as np
import matplotlib.pyplot as plt
import math


def read_data(file_name):
    """Чтение из файла массива"""
    with open(file_name) as inp:
        lines = inp.readlines()
        
    values = []
    for line in lines:
        values.append((float(line.split(" ")[0]), float(line.split(" ")[1])))

    return values


def calc_function(x, values):
    """Вычисление функции для заданного аргумента , где
        x -- аргумент для которого должна быть расчитанна функция
        values -- список значений функции для неизвестных аргументов
    """

    size = len(values)  # размер массива
    f_x = 0  # начальное значение функции

    for k in range(size):
        ckx_no_lambda, lambda_k = 1, 1

        for j in range(size):
            if j != k:
                ckx_no_lambda *= x - values[j][0]  # вычисление Ck(x) без lambda_k
                lambda_k *= values[k][0] - values[j][0]  # вычисление произведения для lambda_k

        f_x += (ckx_no_lambda / lambda_k) * values[k][1]  # вычисление Ln(x)

    return f_x  # значение функции для аргумента x


def printer(matrix):
    """Красивый вывод матрицы"""
    for i in range(len(matrix)):
        print("\t".join([str(round(k, 10)) for k in matrix[i]]))
    print("")


def draw(values):
    """Отрисовка графика"""
    plt.figure(1)
    plt.subplot(211)
    t1 = np.arange(-4.0, 4.0, 0.02)
    func = eval('lambda x: math.exp(-x * x)')

    plt.grid(True)
    # график итерполяционного многочлена
    plt.plot([i[0] for i in values], [i[1] for i in values], 'bo', t1, [calc_function(i, values) for i in t1], 'k', label="atan")
    # график функции
    plt.plot(t1, [func(x) for x in t1], 'k', color="g", label="-x*x")
    plt.legend()
    plt.show()


def main():
    """Основная функция"""

    # чтение и вывод
    values = read_data("input.dat")
    print("Таблица значений функции:")
    printer(values)

    # вычисление приближенного значения (пункт а)
    print("\nВведите точку для определения приближенного значении функции: ")
    point = float(input())
    print("\nПриближенное значение к точке:", )
    print(calc_function(point, values))

    # отрисовка графика y=f(x) и интерполяционного многочлена (пункт б)
    draw(values)


if __name__ == '__main__':
    main()
