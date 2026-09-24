from fractions import Fraction

from .frac import numstr
from .matutil import submat, transpose
from .det import det_minor


def minors_cofactors(M, w):
    n = len(M)
    w.h(f"Миноры и алгебраические дополнения матрицы {n}×{n}", 1)
    w.p("Минор M_ij — определитель матрицы, полученной вычёркиванием i-й строки и j-го столбца. "
        "Алгебраическое дополнение A_ij = (−1)^(i+j)·M_ij.")
    cof = []
    for i in range(n):
        row = []
        for j in range(n):
            w.h(f"Минор M{i + 1}{j + 1}: вычёркиваем строку {i + 1} и столбец {j + 1}", 2)
            mn = submat(M, i, j)
            w.m(f"M{i + 1}{j + 1}", mn)
            dv = det_minor(mn, w, f"M{i + 1}{j + 1}")
            sgn = Fraction(1) if (i + j) % 2 == 0 else Fraction(-1)
            aij = sgn * dv
            row.append(aij)
            w.p(f"A{i + 1}{j + 1} = (−1)^({i + 1}+{j + 1})·M{i + 1}{j + 1} = {numstr(sgn)}·{numstr(dv)} = {numstr(aij)}")
        cof.append(row)

    w.h("Матрица алгебраических дополнений C", 1)
    w.m("C", cof, caption="на каждой позиции (i, j) стоит A_ij", ans=True)

    ad = transpose(cof)
    w.h("Союзная (присоединённая) матрица adj(A) = Cᵀ", 1)
    w.m("adj(A)", ad, caption="получается транспонированием матрицы C", ans=True)
    return cof