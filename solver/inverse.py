from fractions import Fraction

from .frac import latex
from .matutil import transpose, submat
from .det import laplace, det_minor


def _aug_rows(aug):
    n = len(aug[0]) // 2
    rows = []
    for i in range(len(aug)):
        rows.append(aug[i][:n] + [{"sep": True}] + aug[i][n:])
    return rows


def inverse_adjugate(M, w):
    n = len(M)
    w.h("Обратная матрица методом алгебраических дополнений", 1)
    w.m("A", M, caption=f"A — матрица {n}×{n}")

    d = laplace(M, w)
    if d == 0:
        w.p("Так как det(A) = 0, матрица вырожденная — обратной матрицы не существует.", ans=True)
        return None
    w.l(f"\\det(A) = {latex(d)} \\neq 0 \\;\\Rightarrow\\; \\text{{матрица невырожденная}}", ans=True)

    w.h("1) Находим алгебраические дополнения всех элементов", 1)
    cof = []
    for i in range(n):
        row = []
        for j in range(n):
            mn = submat(M, i, j)
            dv = det_minor(mn, w, f"M_{{{i + 1}{j + 1}}}")
            sgn = Fraction(1) if (i + j) % 2 == 0 else Fraction(-1)
            aij = sgn * dv
            row.append(aij)
            w.l(f"A_{{{i + 1}{j + 1}}} = {latex(sgn)} \\cdot M_{{{i + 1}{j + 1}}} = {latex(aij)}")
        cof.append(row)

    w.h("Матрица алгебраических дополнений C", 1)
    w.m("C", cof)
    ad = transpose(cof)
    w.h("2) Союзная матрица adj(A) = Cᵀ", 1)
    w.m("adj(A)", ad)

    w.h("3) A⁻¹ = (1/det(A))·adj(A)", 1)
    w.p(f"Каждый элемент adj(A) делим на det(A) = {latex(d)}:")
    inv = []
    for i in range(n):
        row = []
        for j in range(n):
            v = ad[i][j] / d
            row.append(v)
            w.l(f"(A^{{-1}})_{{{i + 1}{j + 1}}} = \\frac{{{latex(ad[i][j])}}}{{{latex(d)}}} = {latex(v)}")
        inv.append(row)
    w.m("A⁻¹", inv, caption="Обратная матрица", ans=True)
    return inv


def inverse_gauss(M, w):
    n = len(M)
    w.h("Обратная матрица методом Гаусса—Жордана", 1)
    w.m("A", M, caption=f"A — матрица {n}×{n}")
    w.p("Составляем расширенную матрицу [A | E], где E — единичная матрица. "
        "Элементарными преобразованиями строк приводим её к виду [E | A⁻¹].")
    aug = [M[i][:] + [Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]
    w.m("[A | E]", _aug_rows(aug))

    for k in range(n):
        piv = aug[k][k]
        if piv == 0:
            p = None
            for r in range(k + 1, n):
                if aug[r][k] != 0:
                    p = r
                    break
            if p is None:
                raise ValueError("Ведущий элемент не найден — матрица вырождена, обратной не существует.")
            aug[k], aug[p] = aug[p], aug[k]
            w.p(f"Меняем местами строки {k + 1} и {p + 1}.")
            w.m("[A | E]", _aug_rows(aug))
            piv = aug[k][k]

        w.p(f"Нормируем строку {k + 1}: делим на ведущий элемент {latex(piv)}.")
        aug[k] = [v / piv for v in aug[k]]
        w.m("[A | E]", _aug_rows(aug))

        for i in range(n):
            if i == k or aug[i][k] == 0:
                continue
            m = aug[i][k]
            w.l(f"R_{{{i + 1}}} \\leftarrow R_{{{i + 1}}} - \\left({latex(m)}\\right) R_{{{k + 1}}}")
            aug[i] = [aug[i][j] - m * aug[k][j] for j in range(2 * n)]
            w.m("[A | E]", _aug_rows(aug))

    inv = [row[n:] for row in aug]
    w.p("Слева получилась единичная матрица — значит, справа стоит искомая обратная матрица:")
    w.m("A⁻¹", inv, caption="Обратная матрица", ans=True)
    return inv