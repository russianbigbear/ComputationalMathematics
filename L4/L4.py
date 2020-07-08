from numpy import *
import numpy as np

epsilon = 0.0000000001
max_iterations = 10000


def read_data(file_name):
    """Чтение данных их файла"""

    # открытие файла
    with open(file_name) as ifs:
        lines = ifs.readlines()

    # размер матрицы
    size_m = len(lines[0].split())

    # создание пустой матрицы
    a = empty((size_m, size_m), dtype=float_)

    # создание пустого вектора
    vector = empty((size_m, 1), dtype=float_)

    # считывание и формирование  матрицы
    i = 0
    for i in range(size_m):
        lines[i] = lines[i].split()
        for j in range(size_m):
            a[i, j] = float(lines[i][j])

    # считывание и формирование вектора
    lines[i + 2] = lines[i + 2].split()
    for j in range(size_m):
        vector[j, 0] = float(lines[i + 2][j])
    lam0 = float(lines[i + 2][j + 1])

    return a, vector, size_m, lam0


def swap_row(matrix_a, first, second, size_m):
    """Перестановка двух строк"""
    for i in range(size_m):
        x = matrix_a[first][i]
        matrix_a[first][i] = matrix_a[second][i]
        matrix_a[second][i] = x

    return matrix_a


def swap_colomn(matrix_a, first, second, size):
    """Перестановка двух столбцов"""
    for i in range(size):
        x = matrix_a[i][first]
        matrix_a[i][first] = matrix_a[i][second]
        matrix_a[i][second] = x

    return matrix_a


def norma(matrix_a):
    """Вычисление нормы для матрицы"""
    return sqrt(sum(matrix_a * matrix_a))


def find_eigenvalues(a, vector, size_m):
    """Вычисление 2ух наибольших по модулю собственных числа"""

    # Вычисление собственного числа лямда 1 для A
    x = vector[:].transpose().reshape(size_m)
    xn = empty(size_m)  # последовательность, сходящаяся к искому числу

    lam = 0
    last = 1

    i = 0
    # выполняем, пока лямбда не перестанут менятся в пределах заданной точности или в пределах итераций
    while abs(last - lam) >= epsilon and i < max_iterations:
        last = lam
        xn = matmul(a, x)  # произведение координат
        # скалярное произведение, как суииа П соответствующих координат
        lam = sum(xn * x) / sum(x * x)
        x = xn[:]
        i += 1

    iterations_1 = i

    # вычисление собственного вектора для лямда 1
    e_vec_1 = xn / norma(xn)

    # Вычисление собственного числа лямда 2 для At
    at = a.transpose()

    x = vector[:].transpose().reshape(size_m)

    lam2 = 0
    last = 1

    i = 0
    # выполняем, пока лямбда не перестанут менятся в пределах заданной точности или в пределах итераций
    while abs(last - lam2) >= epsilon and i < max_iterations:
        last = lam2
        xn = matmul(at, x)  # произведение координат
        # скалярное произведение, как суииа П соответствующих координат
        lam2 = sum(xn * x) / sum(x * x)
        x = xn[:]
        i += 1

    # вычисление собственного вектора для лямда 2
    g_vec = xn / norma(xn)

    # Вычисление собственного числа лямда 2 для A
    # начальный вектор без пропорциональной компоненты
    y = vector.transpose().reshape(size_m) - (
            (sum(vector.transpose().reshape(size_m) * g_vec)) / (sum(e_vec_1 * g_vec))) * e_vec_1
    yn = empty(size_m)

    lam3 = 0
    last = 1

    i = 0
    # выполняем, пока лямбда не перестанут менятся в пределах заданной точности или в пределах итераций
    while abs(last - lam3) >= epsilon and i < max_iterations:
        last = lam3
        yn = matmul(a, y)  # произведение координат
        # скалярное произведение, как суииа П соответствующих координат
        lam3 = sum(yn * y) / sum(y * y)
        y = yn[:]
        if i % 5 == 0:
            y = y - (sum(y * g_vec)) / (sum(e_vec_1 * g_vec)) * e_vec_1
        i += 1

    # вычисление собственного вектора для лямда 3
    e_vec_2 = yn / norma(yn)

    iterations_2 = i

    #  print("Произведение собственных векторов e_vec_1 и e_vec_2 = %f" % sum(e_vec_1 * e_vec_2))

    return lam, e_vec_1, lam3, e_vec_2, iterations_1, iterations_2


def close_computation(matrix_a, vector, size_m, initial):
    """Вычисление близкого числа к лямбда0(введеному)"""
    b = linalg.inv(matrix_a - initial * identity(size_m))

    x = vector[:].transpose().reshape(size_m)
    lam, com_lambda, last, i = 0, 0, 1, 0

    # считаем обратные итерации со сдвигом
    while abs(last - lam) >= epsilon and i < max_iterations:
        last = lam
        xn = matmul(b, x)
        lam = sum(xn * x) / sum(x * x)
        com_lambda = initial + sum(x * x) / sum(xn * x)
        x = xn[:]
        i += 1

    return com_lambda


def calculating_minimum_eigenvalue(matrix_a, vector, size_m):
    """Вычисление минимального собственного числа для матрицы"""
    c = -9
    b = matrix_a + c * identity(size_m)

    x = vector[:].transpose().reshape(size_m)
    beigen_value = 0
    last = 1
    i = 0

    while abs(last - beigen_value) >= epsilon and i < max_iterations:
        last = beigen_value
        xn = matmul(b, x)
        beigen_value = sum(xn * x) / sum(x * x)
        x = xn[:]
        i += 1

    return beigen_value - c


def main():
    a, vector, size_m, lambda_zero = read_data("input.dat")
    l1, e1, l2, e2, it1, it2 = find_eigenvalues(a, vector, size_m)

    # вывод информации
    print("Первое наибольшее по модулю собственное число А: %f" % l1)
    print("Найдено за %d итераций." % it1)
    print("Его собственный вектор: ", e1.transpose().reshape(size_m))
    print()

    print("Первое наибольшее по модулю собственное число А: %f" % l2)
    print("Найдено за %d итераций." % it2)
    print("Его собственный вектор: ", e2.transpose().reshape(size_m))
    print()

    closest_number = close_computation(a, vector, size_m, lambda_zero)
    print("Собственное число, ближайшее к числу %f, равно %f" % (lambda_zero, closest_number))

    min_eigen_value = calculating_minimum_eigenvalue(a, vector, size_m)
    print("Наименьшее по модулю собственное число матрицы A: %f" % min_eigen_value)


if __name__ == '__main__':
    main()
