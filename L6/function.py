import math


def e(amount):
    """Построение таблицы для экспоненты"""
    step = 8 / amount
    x = - 4
    values = []
    for i in range(amount + 1):
        values.append((x, math.exp(-x * x)))
        x += step
    return values


def atan(amount):
    """Построение таблицы для арктангенса"""
    step = 8 / amount
    x = - 4
    values = []
    for i in range(amount + 1):
        values.append((x, math.atan(50 * x)))
        x += step
    return values


def xxx(amount):
    """Построение таблицы для x^3"""
    step = 8 / amount
    x = - 4
    values = []
    for i in range(amount + 1):
        values.append((x, x * x * x))
        x += step
    return values


def main():
    """Основная функция"""
    with open("input.dat", 'w') as file:
        f_x = atan(5)

        for value in f_x:
            file.write("%f %f\n" % (value[0], value[1]))


if __name__ == '__main__':
    main()
