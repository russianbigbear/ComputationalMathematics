from numpy import dot, identity, sign
from math import sqrt

# Константы
MAX_PRECISION = 0.000001  # точность
MAX_ITER = 1000000  # грань итераций


def is_symmetric(matrix, n):
    """Проверка матрицы на симетричность"""
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True


def get_ukl(k, l, alp, bet, size):
    """Получение Ukl матрицы для поворота"""
    result = identity(size)
    result[k, k] = result[l, l] = alp  # alpha стоит на (k,k) / (l,l) месте
    result[l, k] = bet  # beta стоит на (l,k) месте
    result[k, l] = -bet  # -beta стоит на (k,l) месте
    return result


def printer(matrix):
    """Красивый вывод матрицы"""
    for i in range(len(matrix)):
        print("\t".join([str(round(k, 10)) for k in matrix[i]]))
    print("")


def compare(matrix):
    """Проверка на точность"""
    for l, i in enumerate(matrix):
        for k, j in enumerate(i):
            if k == l:
                continue
            elif abs(j) > MAX_PRECISION:
                return False
    return True


def print_vec(vector, round_value, delim="\t"):
    """Вывод вектора"""
    print(delim.join([str(round(i, round_value)) for i in vector]))


def get_R(matrix, l, e):
    """Вычисление невязки"""
    return [i - j for i, j in zip(dot(matrix, e), [k * l for k in e])]  # rk = bk - sum(aki*xi)


def main():
    """Основная функция"""

    # чтение данных
    with open("input.dat") as file:
        size_A = int(file.readline())  # размер матрицы
        matrix_A = [[float(i) for i in file.readline().split()] for _ in range(size_A)]  # чтение матрицы

    # проверка матрицы
    if not is_symmetric(matrix_A, size_A):
        print("Матрица не симетрична! Исправте входные данные")
        return -1

    copy_A = matrix_A.copy()  # копирование матрицы
    matrix_E = identity(size_A)  # единичная матрица
    count_iterations = 0  # количество итераций
    num_max_string = 0  # номер  максимальной строки
    num_max_element = 1  # номер  максимального элемента

    # вывод матрицы
    print("Введеная матрица:")
    printer(matrix_A)

    # решение
    while count_iterations < MAX_ITER:
        maximum = matrix_A[0][1]

        # поиск элемента для анулирования(максималнього)
        for i, x in enumerate(matrix_A):
            for j, y in enumerate(x):
                if i != j and abs(y) > abs(maximum):
                    maximum, num_max_string, num_max_element = y, i, j

        # если akk = all, то alpha= beta = корень(1/2)
        if abs(matrix_A[num_max_string][num_max_string] - matrix_A[num_max_element][num_max_element]) < MAX_PRECISION:
            alpha = beta = sqrt(0.5)
        # иначе сложное вычисление через mi(u)
        else:
            mi = 2 * matrix_A[num_max_string][num_max_element] / (matrix_A[num_max_string][num_max_string] - matrix_A[num_max_element][num_max_element])
            fraction = 1 / sqrt(1 + mi ** 2)  # дробь в alpha и beta
            alpha = sqrt(0.5 * (1 + fraction))  # alpha
            beta = sign(mi) * sqrt(0.5 * (1 - fraction))  # beta

        matrix_D = get_ukl(num_max_string, num_max_element, alpha, beta, size_A)
        transpose_D = get_ukl(num_max_string, num_max_element, alpha, beta, size_A)
        # транспонирование
        transpose_D[num_max_string][num_max_element], transpose_D[num_max_element][num_max_string] = beta, -beta

        # вывод инфо итерации
        print(count_iterations + 1 , "итерация:")
        # B = U(T)kl * C; C = A * Ukl -> B = U(T)kl * A * Ukl
        B = dot(transpose_D, matrix_A)   # B =  U(T}kl * A
        B = dot(B, matrix_D)  # B = B * Ukl
        print("Матрица B:")
        printer(B)  # вывод

        if compare(B):
            break

        matrix_A = B.copy()  # матрица A = B (вращение)
        matrix_E = dot(matrix_E, matrix_D)  # вычисление собственных векторов

        count_iterations += 1

    # вывод ответа
    print("Ответ: ")
    [print("Собственное число %d: " % (i + 1), matrix_A[i][i]) for i in range(size_A)]
    print("\nКоличество итераций: ", count_iterations + 1)

    print("\nСобственные вектора: \n")

    # вывод векторов
    for i in range(size_A):
        print("Вектор %d:" % (i + 1))
        e = [matrix_E[j][i] for j in range(size_A)]
        print_vec(e, 10)

        print("Невязка:")
        print_vec(get_R(copy_A, matrix_A[i][i], e), 10)
        print("")


if __name__ == '__main__':
    main()
