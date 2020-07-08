import numpy as np
from numpy import *
import copy

epsilon = 0.00001


def read_A(file_name):
    """Чтение файла матрицы"""
    with open(file_name) as ifs:  # Чтение строк
        lines = ifs.readlines()

    m = len(lines[0].split())  # Рабиваем строки
    A = empty((m, m), dtype=float)  # Создание временнной матрицы

    # Создание матрицы А
    for i in range(m):
        lines[i] = lines[i].split()
        for j in range(m):
            A[i, j] = float(lines[i][j])
    return A, m  # Возращаем матрицу A и рамер m


def read_B(file_name, m):
    """Чтение файла ответов"""
    with open(file_name) as ifs:  # Чтение строк
        lines = ifs.readlines()

    b = empty((m, 1), dtype=float)  # Создание временнной матрицы

    # Создание матрицы b(матрицы ответов)
    for i in range(m):
        lines[i] = lines[i].split()
        b[i, 0] = float(lines[i][0])
    return b  # Возращаем матрицу ответов


def swap_row(matrix, first, second, size):
    """Перестановка двух строк"""
    for i in range(size):
        x = matrix[first][i]
        matrix[first][i] = matrix[second][i]
        matrix[second][i] = x
    return matrix


def swap_colomn(matrix, first, second, size):
    """Перестановка двух столбцов"""
    for i in range(size):
        x = matrix[i][first]
        matrix[i][first] = matrix[i][second]
        matrix[i][second] = x
    return matrix


def sd(a, b, size):
    """Вычисление матриц S и D"""
    s = zeros((size, size), dtype=float_)  # Зануление матрицы S
    d = zeros((size, size), dtype=float_)  # Зануление матрицы D
    j = 1
    for i in range(size):
        d[i, i] = np.sign(a[i, i] - (sum((s[:i, i] ** 2) * (d[:i, :i]))))  # Вычисление элемента диагональной  матрицы D
        s[i, i] = np.sqrt(abs(a[i, i] - (sum((s[:i, i] ** 2) * (d[:i, :i])))))  # Вычисление диагонального элемента S
        if abs(s[i, i]) < epsilon:
            swap_row(a, i, j, size)  # Меняем строки в  матрице A
            swap_row(b, i, j, 1)  # Меняем строки в матрице b
            swap_colomn(a, i, j, size)  # Меняем столбцы в A
            j += 1
            i -= 1
            continue
        else:
            j = i

        for k in range(i + 1, size, +1):  # Находим Sik элемент матрицы S
            s[i, k] = (a[i, k] - (sum((s[:i, i]) * (s[:i, k]) * (d[:i, :i])))) / (s[i, i] * d[i, i])

    return s, d  # Возврат полученных матриц


def solve(a, b, size):
    """Вычисление ответа матрицы"""
    s, d = sd(a, b, size)  # Вычисление S и D
    st = np.matrix.transpose(s)  # Вычисление St(трансонирование)
    StDS = np.matmul(np.matmul(np.matrix.transpose(s), np.matrix(d)), np.matrix(s))  # Проверка
    print("Input matrix - A: ")
    print(a)
    print("\nInput answer matrix - B: ")
    print(b)
    print("\nTransposed matrix - St: ")
    print(st)
    print("\nDiagonal matrix - D: ")
    print(d)
    print("\nUpper triangular matrix - S: ")
    print(s)
    print("\nPartition result - StDS: ")
    print(StDS)

    # создание матриц ответов z,y,x
    z = zeros((size, 1))
    y = zeros((size, 1))
    x = zeros((size, 1))

    z[0, 0] = b[0, 0] / s[0, 0]  # z1 = b1/s11 - начальное значение z

    # Нахождение zi
    for i in range(1, size, +1):
        z[i, 0] = (b[i, 0] - sum(z[:i, 0] * s[:i, i])) / s[i, i]

    # Нахождение yi
    for i in range(size):
        y[i, 0] = z[i, 0] / d[i, i]

    x[size - 1, 0] = y[size - 1, 0] / s[size - 1, size - 1]  # xm = ym/smm - начальное значение x

    # Нахождение ответов x
    for i in range(size - 2, -1, -1):
        x[i, 0] = (y[i, 0] - sum(s[i, i + 1:size] * x[i + 1:size, 0])) / s[i, i]
    return x  # Возврат ответов


def norma(matrix, size):
    """Вычисление нормы матрицы"""
    values = []
    for i in range(size):
        values.append(0)
        for j in range(size):
            values[i] += abs(matrix[i, j])
    return max(values)


def gilbert(n):
    """"Создание матрицы Гилюерта и подсчет её числа"""
    h = empty((n, n))
    for i in range(n):
        for j in range(n):
            h[i, j] = 1 / (i + j + 1)

    print("\nGilbert matrix: ")
    print(h)

    return norma(h, n) * norma(linalg.inv(h), n)  # ||A||*||-A||

def main():
    """Основная функция"""
    # Чтение файлов
    a, size = read_A("inputA.dat")
    b = read_B("inputB.dat", size)

    # Копирование матриц
    A_copy = copy.deepcopy(a)
    b_copy = copy.deepcopy(b)

    # Вычисление ответов системы уравнений
    x = solve(a, b, size)  # Решение СЛУ

    if x is False:
        print("The matrix is degenerate!")
    else:
        inv = np.linalg.inv(A_copy)  # Вычисление обратной матрицы
        y = zeros(size)

        # Вычисление невязки
        for i in range(size):
            for j in range(size):
                y[i] += A_copy[i, j] * x[j, 0]
            y[i] = b_copy[i, 0] - y[i]
            # Если величина невязки незначительна, то нулим
            if abs(y[i]) < epsilon:
                y[i] = 0

        # Запись ответов
        with open('Result.dat', 'w') as out:
            out.write("Решение: \n")
            for i in range(size):
                out.write(str(x[i]) + ' ')

            out.write("\nОпределитель\n" + str(np.linalg.det(a)))

            out.write("\nНевязки:\n")
            for i in range(size):
                out.write(str(y[i]) + ' ')

            out.write("\nОбратная матрица: \n")
            for i in range(size):
                for j in range(size):
                    out.write(str(inv[i, j]) + ' ')
                out.write('\n')

        out.close()

if __name__ == '__main__':
    main()
