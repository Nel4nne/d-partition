"""
Простая итерационная программа для построения визуализации
D-разбиения для уравнения eps*x'(t) + x(t) + a1*x(t-p) + a2*x(t-q) = 0
"""

from matplotlib import pyplot
from math import sin, cos


def calc_a1_wrapper(p: float, q: float, eps: float) -> callable:
    def calc_a1(w: float) -> float:
        try:
            return (sin(w*q) - eps*w*cos(w*q)) / sin(w*(p-q))
        except ZeroDivisionError as zde:
            print(f"Возникла ошибка: деление на 0 при значении w = {w}")
            return None

    return calc_a1


def calc_a2_wrapper(p: float, q: float, eps: float) -> callable:
    def calc_a2(w: float) -> float:
        try:
            return (eps*w*cos(w*p) - sin(w*p)) / sin(w*(p-q))
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
        func_points.append(func(arg))
        args.append(arg)
        arg += arg_step

    return args, func_points


def main():
    # Задаём параметры по условию 
    p = 1
    q = 5
    eps = 0.02
    w_start = 0
    w_end = 100000
    precision = 1e6

    figure, (ax1, ax2) = pyplot.subplots(nrows=1, ncols=2)

    # Инициализируем функции параметрами из условия
    calc_a1 = calc_a1_wrapper(p, q, eps)
    calc_a2 = calc_a2_wrapper(p, q, eps)

    # Табулируем функции по 'w'
    w1_list, a1_list = tab_func(calc_a1, w_start, w_end, precision)
    w2_list, a2_list = tab_func(calc_a2, w_start, w_end, precision)
    
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