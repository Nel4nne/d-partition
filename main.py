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
from math import sin, cos, pi


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


def _calc_lambda_1_wrapper(p: float, q: float, eps: float) -> callable:
    def _calc_lambda_1(a1: float, a2: float, w: float):
        return a1 * p + a2 * q * cos(w*(p-q)) - eps * cos(w*p)
    
    return _calc_lambda_1


def _calc_lambda_2_wrapper(p: float, q: float, eps: float) -> callable:
    def _calc_lambda_2(a1: float, a2: float, w: float):
        return a2 * q + a1 * p * cos(w*(p-q)) - eps * cos(w*q)
    
    return _calc_lambda_2


def _calc_eps_zero_wrapper(p: float, q: float, w: float):
    def _calc_eps_zero(a1: float):
        ep = 1e-6
        return (-1 - a1 * cos(w*p))/cos(w*q) if abs(cos(w*q)) > ep else nan

    return _calc_eps_zero


def _calc_lambda_zero(a1: float):
    return -a1 - 1


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


def is_zero(val: float) -> bool:
    eps = 1e-10
    return abs(val) < eps


def _calc_root_transition(
    calc_l1_func: callable,
    calc_l2_func: callable,
    w: list,
    a1: list,
    a2: list
) -> tuple[list, list]:
    eps_zero = 1e-10 # для сравнения лямбды с нулём
    a1_shift = 0.2 # для отступа от кривой
    transition_a1_list = []
    transition_a2_list = []
    # Для оставшихся точек считаем знак lambda_1 и lambda_2 (если потребуется)
    # чтобы понять в какую комплексную полуплоскость переходит корень
    for point in zip(w, a1, a2):
        # Выколотые пропускаем
        if any((point[1] is None, point[2] is None)):
            transition_a1_list.append(nan)
            transition_a2_list.append(nan)
            continue
        
        # Считаем лямбду в точке на кривой
        lambda_root = calc_l1_func(point[1], point[2], point[0])
        # Если равна нулю то нужно считать вторую
        if is_zero(lambda_root):
            lambda_root = calc_l2_func(point[1], point[2], point[0])
            # Если и второй ноль, то дальнейшие приближения пусть считает тот, кому надоело жить
            # а я это место просто пропущу...
            if is_zero(lambda_root):
                transition_a1_list.append(nan)
                transition_a2_list.append(nan)
                continue
            # Если отрицательное, то корень переходит в левую комплексную полуплоскость,
            # ставим точку справа от кривой
            if lambda_root < 0:
                transition_a2_list.append(point[2] + a1_shift)
            # А если положительное лямбда, то ставим точку слева от кривой
            else:
                transition_a2_list.append(point[2] - a1_shift)
            transition_a1_list.append(point[1])
        else:
            # Если отрицательное, то корень переходит в левую комплексную полуплоскость,
            # ставим точку справа от кривой
            if lambda_root < 0:
                transition_a1_list.append(point[1] + a1_shift)
            # А если положительное лямбда, то ставим точку слева от кривой
            else:
                transition_a1_list.append(point[1] - a1_shift)
            transition_a2_list.append(point[2])
    
    return transition_a1_list, transition_a2_list


def main():
    # Задаём параметры по условию 
    p = 1
    q = 5
    eps = 0.5
    w_start = 0.0001
    w_end = 20
    precision = 1e6

    # Инициализируем функции параметрами из условия
    calc_a1 = calc_a1_wrapper(p, q, eps)
    calc_a2 = calc_a2_wrapper(p, q, eps)
    # Инициализация функций для подсчета перехода корней
    calc_lambda_1 = _calc_lambda_1_wrapper(p, q, eps)
    calc_lambda_2 = _calc_lambda_2_wrapper(p, q, eps)

    # Расчёт по основному решению (когда определитель не 0)
    # Табулируем функции по 'w'
    w1_list, a1_list = tab_func(calc_a1, w_start, w_end, precision)
    w2_list, a2_list = tab_func(calc_a2, w_start, w_end, precision)
    # Удаляем все точки, которые не попадают в нужный прямоугольник
    exclude_inf_points(
        w1_list,
        w2_list,
        a1_list,
        a2_list,
        (-10, 10, 10, -10) # Левый - правый - верхний - нижний
    )
    # Для отображения перехода корней для основного решения
    (
        transition_roots_a1,
        transition_roots_a2
    ) = _calc_root_transition(
        calc_lambda_1,
        calc_lambda_2,
        w1_list,
        a1_list,
        a2_list,
    )

    # Расчёт для нулевого решения
    # Нахождение точек (a1, a2) для нулевого решения
    a1_start = -10
    a1_end = 10
    lambda_zero_a1, lambda_zero_a2 = tab_func(_calc_lambda_zero, a1_start, a1_end, precision)

    figure, ax2 = pyplot.subplots(nrows=1, ncols=1)
    # Расчёт для решений, когда определитель системы равен нулю
    # Если определитель системы равен нулю
    for k in range(5):
        w = (pi * k)/(q - p)
        if abs(cos(w*q)) < 1e-6:
            continue
        calc_eps_zero = _calc_eps_zero_wrapper(p, q, w)
        eps_zero_a1, eps_zero_a2 = tab_func(calc_eps_zero, a1_start, a1_end, precision)
        (
            transition_roots_det_zero_a1,
            transition_roots_det_zero_a2,
        ) = _calc_root_transition(
            calc_lambda_1,
            calc_lambda_2,
            [w] * len(eps_zero_a1),
            eps_zero_a1,
            eps_zero_a2
        )
        ax2.plot(eps_zero_a1, eps_zero_a2, color='b')
        # ax2.plot(transition_roots_det_zero_a1, transition_roots_det_zero_a2, color='g')

    # Строим общий график(портрет)
    ax2.set_xlabel('a1')
    ax2.set_ylabel('a2')
    ax2.plot(a1_list, a2_list)
    # ax2.plot(transition_roots_a1, transition_roots_a2, color='g')
    ax2.plot(lambda_zero_a1, lambda_zero_a2)
    # ax2.plot(eps_zero_a1, eps_zero_a2)

    pyplot.show()


if __name__ == "__main__":
    main()