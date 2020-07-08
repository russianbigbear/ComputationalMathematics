from numpy import *
import numpy as np

epsilon = 0.00001


def read_data(file_name):
    """Чтение даных из файла"""
    with open(file_name) as ifs:  # Чтение строк
        lines = ifs.readlines()

    m = len(lines[0].split()) - 1  # Рабиваем строки

    a = empty((m, m), dtype=float_)  # Создание временной матрицы A
    b = empty((m, 1), dtype=float_)  # Создание временной матрицы b

    # Заполнение матрицы А
    for i in range(m):
        lines[i] = lines[i].split()
        for j in range(m):
            a[i, j] = float(lines[i][j])

        b[i, 0] = float(lines[i][m])  # Заполнение матрицы b
    return a, b, m  # Возвращаем данные


def swap_row(matrix, first, second, size):
    """Перестановка строк"""
    for i in range(size):
        x = matrix[first][i]
        matrix[first][i] = matrix[second][i]
        matrix[second][i] = x
    return matrix


def swap_colomn(matrix, first, second, size):
    """Перестановка столбцов"""
    for i in range(size):
        x = matrix[i][first]
        matrix[i][first] = matrix[i][second]
        matrix[i][second] = x
    return matrix


def norma(matrix, size):
    """Вычисление нормы B"""
    values = []
    for i in range(size):
        values.append(0)
        for j in range(size):
            values[i] += abs(matrix[i, j])
    return max(values)


def fix(a, size):
    """Избавление от диагональных нулей/проверка"""
    for i in range(size):
        if a[i, i] == 0:
            for j in range(i + 1, size, +1):
                if a[i, j] != 0:
                    swap_row(a, i, j, size)
                    break
            if a[i, i] == 0:
                return False
    return True


def convergence(a, size):
    """Проверка дсотаточного условия сходимости метода Якоби"""
    for i in range(size):
        s = 0
        for j in range(size):
            if i != j:
                s += abs(a[i, j])
        if abs(a[i, i]) <= s:
            return False
    return True



def solve_yakobi(a, b, size):
    """Решение системы линейных уравнений"""
    x_new = zeros(size)
    x_old = zeros(size)
    iterations = 0

    while True:
        iterations += 1
        for i in range(size):
            s = 0
            for j in range(size):
                if i != j:
                    s += (a[i, j] / a[i, i]) * x_old[j]  # Сумма (aij/aii)*xj
            x_new[i] = (b[i, 0] / a[i, i]) - s  # Находим xi

        if max(x_new - x_old) < epsilon:
            print("Количество итераций: ") # Выводим количество итераций
            print(iterations)
            return x_new
        else:
            x_old = x_new[:]  # Иначе продолжим итерационный метод


def solve_zeidel(a, b, size):
    """Решение системы линейных уравнений"""
    x_new = zeros(size)
    x_old = zeros(size)

    conv = False
    iterations = 0
    while not conv:
        x_new = copy(x_old)
        iterations += 1
        for i in range(size):
            s1 = sum(a[i, j] * x_new[j] for j in range(i))  # Считаем сумму до диагонального элемента
            s2 = sum(a[i, j] * x_old[j] for j in range(i + 1, size))  # Считаем сумму после диагонального элемента
            x_new[i] = (b[i] - s1 - s2) / a[i, i]  # Вычитаем полученнное и делим на диагональный

        conv = max(x_new - x_old) <= epsilon
        x_old = copy(x_new)

    print("Количество итераций: ")  # Выводим количество итераций
    print(iterations)
    return x_new


def main():
    """Основная функция"""
    a, b, size = read_data("input.dat")

    #  Копирование A и b
    a_copy1 = a_copy = a1 = copy(a)
    b_copy1 = b_copy = b1 = copy(b)

    if fix(a, size) is False:
        print("На главной диагонали присутствует 0, который не убрать.")
        return 1

    if convergence(a, size) is False:
        print("Достаточное условие сходимости не выполнено")
        return 1

    # Решение системы, если все условия выполнены
    print("Решение методом Якоби: ")
    x = solve_yakobi(a1,  b1, size)
    print(x)

    z = zeros(size)
    for i in range(size):
        for j in range(size):
            z[i] += a_copy1[i, j] * x[j]
        z[i] = b_copy1[i, 0] - z[i]

        if abs(z[i]) < epsilon:
            z[i] = 0

    print("Невязки:")
    print(z)

    print("\n")

    print("Решение методом Зейделя: ")
    y = solve_zeidel(a, b, size)
    print(y)

    z1 = zeros(size)
    for i in range(size):
        for j in range(size):
            z1[i] += a_copy[i, j] * y[j]
        z1[i] = b_copy[i, 0] - z1[i]

        if abs(z1[i]) < epsilon:
            z1[i] = 0

    print("Невязки:")
    print(z1)


if __name__ == '__main__':
    main()