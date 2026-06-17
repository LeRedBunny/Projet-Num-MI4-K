'''Functions to approximate the first and second derivatives of a function'''

from types import Interval, Function



def derivative (function: Function, x_int: Interval) -> Function :
    '''Approximates the derivative of f, returns an empty array if failed.'''
    if len(function) != len(x_int) :
        return []

    result = [(function[1] - function[0]) / (x_int[1] - x_int[0])]
    for i in range(1, len(function) - 1) :
        result.append((function[i + 1] - function[i - 1]) / (x_int[i + 1] - x_int[i - 1]))
    result.append((function[-1] - function[-2]) / (x_int[-1] - x_int[-2]))

    return result


def derivative2 (function: Function, x_int: Interval) -> Function :
    '''Approximates the second derivative of f, returns an empty array if failed.'''
    return derivative(derivative(function, x_int), x_int)



# Tests

if __name__ == '__main__' :

    delta = 1e-1
    x_0 = 0
    n = 100

    tab_x = [x_0 + i * delta for i in range(n)]
    tab_f = [(x_0 + i * delta) ** 2 for i in range(n)]
    
    tab_derivee = [2 * x for x in tab_x[:-1]]
    tab_approximation = derivative(tab_f, tab_x)

    tab_derivee2 = [2 for x in tab_x[:-2]]
    tab_approximation2 = derivative2(tab_f, tab_x)

    tab_erreur = [tab_approximation2[i] - tab_derivee2[i] for i in range(n - 2)]

    print(f'Erreur moyenne : {round(sum(tab_erreur) / n, 3)}')
