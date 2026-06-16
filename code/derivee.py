

def derivee (tab_f: list[float], tab_x: list[float]) -> list[float] :
    '''Calcule la dérivée de f, renvoie un tableau vide si échec'''
    if len(tab_f) != len(tab_x) :
        return []
    return [(tab_f[i + 1] - tab_f[i]) / (tab_x[i + 1] - tab_x[i]) for i in range(len(tab_f) - 1)]


def derivee2 (tab_f: list[float], tab_d: list[float], tab_x: list[float]) -> list[float] :
    '''Calcule la dérivée seconde de f, renvoie un tableau vide si échec'''
    if len(tab_f) != len(tab_x) or len(tab_f) != len(tab_d) + 1 :
        return []
    return [(tab_f[i + 2] - 2 * tab_f[i + 1] + tab_f[i]) / (tab_x[i + 1] - tab_x[i]) ** 2 for i in range(len(tab_f) - 2)]



# Tests
if __name__ == '__main__' :

    delta = 1e-1
    x_0 = 0
    n = 100

    tab_x = [x_0 + i * delta for i in range(n)]
    tab_f = [(x_0 + i * delta) ** 2 for i in range(n)]
    
    tab_derivee = [2 * x for x in tab_x[:-1]]
    tab_approximation = derivee(tab_f, tab_x)

    tab_derivee2 = [2 for x in tab_x[:-2]]
    tab_approximation2 = derivee2(tab_f, tab_approximation, tab_x)

    tab_erreur = [tab_approximation2[i] - tab_derivee2[i] for i in range(n - 2)]

    print(f'Erreur moyenne : {round(sum(tab_erreur) / n, 3)}')
