"""
Простая итерационная программа для построения визуализации
D-разбиения для уравнения eps*x'(t) + x(t) + a1*x(t-p) + a2*x(t-q) = 0

Обозначения:
    - 'w' греческая буква омега, при построении D-разбиеня является параметризующим
    - 'eps' малый параметр, где 0 < eps << 1
    - 'p', 'q' данные из условия коэффициенты запаздывания
"""

from matplotlib import pyplot
from numpy import nan
from math import sin, cos


def calc_a1_wrapper(p: float, q: float, eps: float) -> callable:
    def calc_a1(w: float) -> float|None:
        try:
            # Тут задается функция параметризованная от 'w' для подсчета a1
            return (-sin(w*q) - eps*w*cos(w*q)) / sin(w*(q-p))
        except ZeroDivisionError as zde:
            print(f"Возникла ошибка: деление на 0 при значении w = {w}")
            return None

    return calc_a1


def calc_a2_wrapper(p: float, q: float, eps: float) -> callable:
    def calc_a2(w: float) -> float|None:
        try:
            # Тут задается функция параметризованная от 'w' для подсчета a2
            return (eps*w*cos(w*p) + sin(w*p)) / sin(w*(q-p))
        except ZeroDivisionError as zde:
            print(f"Возникла ошибка: деление на 0 при значении w = {w}")
            return None
    
    return calc_a2


def tab_func(func: callable, arg_start: float, arg_end: float, precision: float) -> list[float]:
    arg_step = (arg_end - arg_start) / precision
    arg = arg_start
    func_points = []
    args = []
    while arg_start <= arg <= arg_end:
        func_result = func(arg)
        func_points.append(func(arg))
        args.append(arg)
        arg += arg_step

    return args, func_points


def _point_in_square(
    point: tuple[float, float],
    square: tuple[float, float, float, float]
) -> bool:
    return (
        (point[0] >= square[0] and point[0] <= square[2]) and # Попадает по первой координате 
        (point[1] <= square[1] and point[1] >= square[3]) # Попадает по 2 координате
    )


def exclude_inf_points(
    w1_list: list[float],
    w2_list: list[float],
    a1_list: list[float],
    a2_list: list[float],
    square: tuple[float, float, float, float],
) -> None:
    i = 0
    while i < len(a1_list):
        # Если какая-то из координат не попадает в прямоугольник, то
        # она исключается из всех списков
        if not _point_in_square((a1_list[i], a2_list[i]), square):
            w1_list[i], w2_list[i], a1_list[i], a2_list[i] = [nan] * 4
        i += 1


def main():
    # Задаём параметры по условию 
    p = 1
    q = 2
    eps = 0.7
    w_start = 0.0001
    w_end = 50
    precision = 1e6

    # Инициализируем функции параметрами из условия
    calc_a1 = calc_a1_wrapper(p, q, eps)
    calc_a2 = calc_a2_wrapper(p, q, eps)

    # Табулируем функции по 'w'
    w1_list, a1_list = tab_func(calc_a1, w_start, w_end, precision)
    w2_list, a2_list = tab_func(calc_a2, w_start, w_end, precision)

    # Удаляем все точки, которые не попадают в нужный прямоугольник
    exclude_inf_points(
        w1_list,
        w2_list,
        a1_list,
        a2_list,
        (-100, 100, 100, -100) # Левый - правый - верхний - нижний
    )

    # Создаю два плота для двух графиков
    figure, (ax1, ax2) = pyplot.subplots(nrows=1, ncols=2)

    # Строим график для фукций a1(w) и a2(w) которые являются границами
    # множеств Dn
    ax1.set_xlabel('w')
    ax1.set_ylabel('a1, a2')
    ax1.plot(w1_list, a1_list, label='a1(w)')
    ax1.plot(w2_list, a2_list, label='a2(w)')
    ax1.legend()

    # Строим общий график(портрет)
    ax2.set_xlabel('a1')
    ax2.set_ylabel('a2')
    ax2.plot(a1_list, a2_list)

    pyplot.show()


if __name__ == "__main__":
    main()