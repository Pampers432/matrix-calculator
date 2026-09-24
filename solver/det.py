from fractions import Fraction

from .frac import numstr
from .matutil import submat


def render_terms(terms):
    parts = []
    total = Fraction(0)
    for i, (sgn, val) in enumerate(terms):
        tv = sgn * val
        total += tv
        if tv >= 0:
            parts.append(("" if i == 0 else " + ") + numstr(tv))
        else:
            parts.append((" − " if i else "−") + numstr(-tv))
    return "".join(parts) if parts else "0", total


def det2_formula(M, w, label="det"):
    a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
    p1 = a * d
    p2 = b * c
    val = p1 - p2
    w.p(f"{label} = a·d − b·c = {numstr(a)}·{numstr(d)} − {numstr(b)}·{numstr(c)} = {numstr(p1)} − {numstr(p2)} = {numstr(val)}")
    return val


def sarrus(M, w, label="det"):
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
    w.p(f"t1 = a11·a22·a33 = {numstr(M00)}·{numstr(M11)}·{numstr(M22)} = {numstr(v1)}")
    w.p(f"t2 = a12·a23·a31 = {numstr(M01)}·{numstr(M12)}·{numstr(M20)} = {numstr(v2)}")
    w.p(f"t3 = a13·a21·a32 = {numstr(M02)}·{numstr(M10)}·{numstr(M21)} = {numstr(v3)}")
    w.p(f"t4 = a13·a22·a31 = {numstr(M02)}·{numstr(M11)}·{numstr(M20)} = {numstr(v4)}")
    w.p(f"t5 = a12·a21·a33 = {numstr(M01)}·{numstr(M10)}·{numstr(M22)} = {numstr(v5)}")
    w.p(f"t6 = a11·a23·a32 = {numstr(M00)}·{numstr(M12)}·{numstr(M21)} = {numstr(v6)}")
    line, total = render_terms([
        (Fraction(1), v1),
        (Fraction(1), v2),
        (Fraction(1), v3),
        (Fraction(-1), v4),
        (Fraction(-1), v5),
        (Fraction(-1), v6),
    ])
    w.p(f"{label} = t1 + t2 + t3 − t4 − t5 − t6 = {line} = {numstr(total)}")
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
        w.p(f"{label} = {numstr(M[0][0])}")
        return M[0][0]
    if n == 2:
        return det2_formula(M, w, label)
    if n == 3:
        return sarrus(M, w, label)
    val = gauss_value(M)
    w.p(f"{label} = {numstr(val)}  (значение, найдено методом Гаусса)")
    return val


def laplace(M, w, expand=("row", 0), label="A", heading=True):
    n = len(M)
    if heading:
        w.h(f"Определитель матрицы {label} (матрица {n}×{n})", 1)
    return _expand(M, w, expand, label)


def _expand(M, w, expand, label):
    n = len(M)
    if n == 1:
        w.p(f"det({label}) = {numstr(M[0][0])}")
        return M[0][0]
    if n == 2:
        return det2_formula(M, w, f"det({label})")

    kind, idx = expand
    if kind == "row":
        w.h(f"Разложение по строке {idx + 1}", 2)
    else:
        w.h(f"Разложение по столбцу {idx + 1}", 2)
    w.p("det(A) = Σ (−1)^(i+j)·a_ij·M_ij, где M_ij — минор (определитель, полученный вычёркиванием строки i и столбца j)")

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
        w.p(f"Элемент a{i + 1}{j + 1} = {numstr(el)}:")
        mn = submat(M, i, j)
        w.m(f"Минор M{i + 1}{j + 1}", mn,
            caption=f"вычёркиваем строку {i + 1} и столбец {j + 1}")
        dv = det_minor(mn, w, f"M{i + 1}{j + 1}")
        sgn = Fraction(-1) ** (i + j + 2)
        cof = sgn * dv
        w.p(f"Алгебраическое дополнение A{i + 1}{j + 1} = (−1)^({i + 1}+{j + 1})·M{i + 1}{j + 1} = {numstr(sgn)}·{numstr(dv)} = {numstr(cof)}")
        w.p(f"Слагаемое a{i + 1}{j + 1}·A{i + 1}{j + 1} = {numstr(el)}·{numstr(cof)} = {numstr(el * cof)}")
        terms.append((sgn, el * dv))

    line, total = render_terms(terms)
    w.p(f"det({label}) = {line} = {numstr(total)}", ans=True)
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
            w.p(f"R{i + 1} ← R{i + 1} − ({numstr(m)})·R{k + 1}")
            A[i] = [A[i][j] - m * A[k][j] for j in range(n)]
            w.m("Матрица после преобразования", A)

    prods = [A[k][k] for k in range(n)]
    prod = Fraction(1)
    for x in prods:
        prod *= x
    det = sign * prod
    w.p("Матрица приведена к треугольному виду. Произведение диагональных элементов: "
        + " · ".join(numstr(x) for x in prods) + f" = {numstr(prod)}")
    if swaps:
        w.p(f"Было {swaps} перестановок строк, поэтому знак меняется:  det = (−1)^{swaps}·произведение = {numstr(det)}", ans=True)
    else:
        w.p(f"Перестановок строк не было, поэтому det = произведение диагональных элементов = {numstr(det)}", ans=True)
    return det