from fractions import Fraction

from .frac import latex, ljoin
from .matutil import submat


def render_terms(terms):
    parts = []
    total = Fraction(0)
    for i, (sgn, val) in enumerate(terms):
        tv = sgn * val
        total += tv
        if tv >= 0:
            parts.append((" + " if i else "") + latex(tv))
        else:
            parts.append((" - " if i else "-") + latex(-tv))
    return "".join(parts) if parts else "0", total


def det2_formula(M, w, label="\\det(A)"):
    a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
    p1 = a * d
    p2 = b * c
    val = p1 - p2
    w.l(f"{label} = a \\cdot d - b \\cdot c = {latex(a)} \\cdot {latex(d)} - {latex(b)} \\cdot {latex(c)} = {latex(p1)} - {latex(p2)} = {latex(val)}")
    return val


def sarrus(M, w, label="\\det"):
    M00, M01, M02 = M[0]
    M10, M11, M12 = M[1]
    M20, M21, M22 = M[2]
    v1 = M00 * M11 * M22
    v2 = M01 * M12 * M20
    v3 = M02 * M10 * M21
    v4 = M02 * M11 * M20
    v5 = M01 * M10 * M22
    v6 = M00 * M12 * M21
    w.p("Правило Саррюса (правило треугольников):")
    w.l(f"t_1 = a_{{11}} a_{{22}} a_{{33}} = {latex(M00)} \\cdot {latex(M11)} \\cdot {latex(M22)} = {latex(v1)}")
    w.l(f"t_2 = a_{{12}} a_{{23}} a_{{31}} = {latex(M01)} \\cdot {latex(M12)} \\cdot {latex(M20)} = {latex(v2)}")
    w.l(f"t_3 = a_{{13}} a_{{21}} a_{{32}} = {latex(M02)} \\cdot {latex(M10)} \\cdot {latex(M21)} = {latex(v3)}")
    w.l(f"t_4 = a_{{13}} a_{{22}} a_{{31}} = {latex(M02)} \\cdot {latex(M11)} \\cdot {latex(M20)} = {latex(v4)}")
    w.l(f"t_5 = a_{{12}} a_{{21}} a_{{33}} = {latex(M01)} \\cdot {latex(M10)} \\cdot {latex(M22)} = {latex(v5)}")
    w.l(f"t_6 = a_{{11}} a_{{23}} a_{{32}} = {latex(M00)} \\cdot {latex(M12)} \\cdot {latex(M21)} = {latex(v6)}")
    line, total = render_terms([
        (Fraction(1), v1),
        (Fraction(1), v2),
        (Fraction(1), v3),
        (Fraction(-1), v4),
        (Fraction(-1), v5),
        (Fraction(-1), v6),
    ])
    w.l(f"{label} = t_1 + t_2 + t_3 - t_4 - t_5 - t_6 = {line} = {latex(total)}")
    return total


def gauss_value(M):
    A = [r[:] for r in M]
    n = len(A)
    sign = Fraction(1)
    for k in range(n):
        piv = A[k][k]
        if piv == 0:
            p = None
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    p = r
                    break
            if p is None:
                return Fraction(0)
            A[k], A[p] = A[p], A[k]
            sign = -sign
            piv = A[k][k]
        for i in range(k + 1, n):
            if A[i][k] != 0:
                m = A[i][k] / piv
                A[i] = [A[i][j] - m * A[k][j] for j in range(n)]
    prod = Fraction(1)
    for k in range(n):
        prod *= A[k][k]
    return sign * prod


def det_minor(M, w, label):
    n = len(M)
    if n == 1:
        w.l(f"{label} = {latex(M[0][0])}")
        return M[0][0]
    if n == 2:
        return det2_formula(M, w, label)
    if n == 3:
        return sarrus(M, w, label)
    val = gauss_value(M)
    w.l(f"{label} = {latex(val)}  \\text{{(значение, найдено методом Гаусса)}}")
    return val


def laplace(M, w, expand=("row", 0), label="A", heading=True):
    n = len(M)
    if heading:
        w.h(f"Определитель матрицы {label} (матрица {n}×{n})", 1)
    return _expand(M, w, expand, label)


def _expand(M, w, expand, label):
    n = len(M)
    if n == 1:
        w.l(f"\\det({label}) = {latex(M[0][0])}")
        return M[0][0]
    if n == 2:
        return det2_formula(M, w, f"\\det({label})")

    kind, idx = expand
    if kind == "row":
        w.h(f"Разложение по строке {idx + 1}", 2)
    else:
        w.h(f"Разложение по столбцу {idx + 1}", 2)
    w.l("\\det(A) = \\sum (-1)^{i+j}\\, a_{ij} \\, M_{ij}")
    w.p("Здесь M_ij — минор (определитель матрицы, полученной вычёркиванием строки i и столбца j).")

    terms = []
    for k in range(n):
        if kind == "row":
            i, j = idx, k
        else:
            i, j = k, idx
        el = M[i][j]
        if el == 0:
            w.p(f"Элемент a{i + 1}{j + 1} = 0 — слагаемое равно нулю, пропускаем.")
            continue
        w.p(f"Элемент a{i + 1}{j + 1} = {latex(el)}:")
        mn = submat(M, i, j)
        w.m(f"M_{{{i + 1}{j + 1}}}", mn,
            caption=f"вычёркиваем строку {i + 1} и столбец {j + 1}")
        dv = det_minor(mn, w, f"M_{{{i + 1}{j + 1}}}")
        sgn = Fraction(-1) ** (i + j + 2)
        cof = sgn * dv
        w.l(f"A_{{{i + 1}{j + 1}}} = (-1)^{{{i + 1 + j + 1}}} \\, M_{{{i + 1}{j + 1}}} = {latex(sgn)} \\cdot \\left({latex(dv)}\\right) = {latex(cof)}")
        w.l(f"a_{{{i + 1}{j + 1}}} A_{{{i + 1}{j + 1}}} = {latex(el)} \\cdot \\left({latex(cof)}\\right) = {latex(el * cof)}")
        terms.append((sgn, el * dv))

    line, total = render_terms(terms)
    w.l(f"\\det({label}) = {line} = {latex(total)}", ans=True)
    return total


def det_gauss(M, w):
    n = len(M)
    w.h(f"Определитель матрицы (матрица {n}×{n}) методом Гаусса", 1)
    A = [r[:] for r in M]
    w.m("A", A, caption=f"Исходная матрица {n}×{n}")
    sign = Fraction(1)
    swaps = 0
    for k in range(n):
        piv = A[k][k]
        if piv == 0:
            p = None
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    p = r
                    break
            if p is None:
                w.p("Ведущий (главный) элемент нулевой и ненулевых элементов под ним нет — определитель равен 0.", ans=True)
                return Fraction(0)
            A[k], A[p] = A[p], A[k]
            sign = -sign
            swaps += 1
            w.p(f"Меняем местами строки {k + 1} и {p + 1} (при перестановке строк определитель меняет знак).")
            w.m("Матрица после перестановки строк", A)
            piv = A[k][k]

        for i in range(k + 1, n):
            if A[i][k] == 0:
                continue
            m = A[i][k] / piv
            w.l(f"R_{{{i + 1}}} \\leftarrow R_{{{i + 1}}} - \\left({latex(m)}\\right) R_{{{k + 1}}}")
            A[i] = [A[i][j] - m * A[k][j] for j in range(n)]
            w.m("Матрица после преобразования", A)

    prods = [A[k][k] for k in range(n)]
    prod = Fraction(1)
    for x in prods:
        prod *= x
    det = sign * prod
    w.l("\\det(A) = " + " \\cdot ".join(latex(x) for x in prods) + f" = {latex(prod)}")
    if swaps:
        w.l(f"\\det(A) = (-1)^{{{swaps}}} \\cdot \\prod a_{{{{ii}}}} = {latex(det)}", ans=True)
    else:
        w.l(f"\\det(A) = \\prod a_{{{{ii}}}} = {latex(det)}", ans=True)
    return det