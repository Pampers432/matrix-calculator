from fractions import Fraction

from .frac import latex
from .matutil import transpose


def _par(v):
    return f"\\left({latex(v)}\\right)" if v < 0 else latex(v)


def add_sub(A, B, w, op="+"):
    name = {"+": "Сложение матриц", "-": "Вычитание матриц"}[op]
    w.h(f"{name}: C = A {op} B", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}")
    w.m("B", B, caption=f"B — матрица {m}×{n}")
    w.p("Размерности совпадают, поэтому складываем (вычитаем) элементы с одинаковыми индексами.")
    C = []
    for i in range(m):
        row = []
        for j in range(n):
            v = A[i][j] + B[i][j] if op == "+" else A[i][j] - B[i][j]
            row.append(v)
            sig = " + " if op == "+" else " - "
            nums = f"{_par(A[i][j])}{sig}{_par(B[i][j])}"
            w.l(f"c_{{{i + 1}{j + 1}}} = a_{{{i + 1}{j + 1}}} "
                f"{'+ ' if op == '+' else '- '}b_{{{i + 1}{j + 1}}} = {nums} = {latex(v)}")
        C.append(row)
    w.m(f"C = A {op} B", C, caption="Результат", ans=True)
    return C


def scalar_mul(k, A, w):
    w.h("Умножение матрицы на число", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}: умножаем каждый элемент на k = {latex(k)}")
    C = []
    for i in range(m):
        row = []
        for j in range(n):
            v = k * A[i][j]
            row.append(v)
            w.l(f"c_{{{i + 1}{j + 1}}} = {latex(k)} \\cdot a_{{{i + 1}{j + 1}}} = {latex(k)} \\cdot \\left({latex(A[i][j])}\\right) = {latex(v)}")
        C.append(row)
    w.m("C = k·A", C, caption="Результат", ans=True)
    return C


def _mt(a, b):
    ta = f"\\left({latex(a)}\\right)" if a < 0 else latex(a)
    tb = f"\\left({latex(b)}\\right)" if b < 0 else latex(b)
    return f"{ta} \\cdot {tb}"


def mat_mul(A, B, w):
    w.h("Умножение матриц: C = A·B", 1)
    m, p = len(A), len(A[0])
    n = len(B[0])
    w.m("A", A, caption=f"A — матрица {m}×{p}")
    w.m("B", B, caption=f"B — матрица {p}×{n}")
    w.p(f"Число столбцов A ({p}) равно числу строк B ({p}) — матрицы согласованы, "
        f"результат C имеет размер {m}×{n}.")
    w.p("Элемент c_ij равен сумме произведений элементов i-й строки A на j-й столбец B.")
    C = []
    for i in range(m):
        row = []
        for j in range(n):
            terms = [A[i][k] * B[k][j] for k in range(p)]
            total = sum(terms, Fraction(0))
            row.append(total)
            formula = " + ".join(f"a_{{{i + 1}{k + 1}}} \\cdot b_{{{k + 1}{j + 1}}}" for k in range(p))
            products = " + ".join(_mt(A[i][k], B[k][j]) for k in range(p))
            w.l(f"c_{{{i + 1}{j + 1}}} = {formula} = {products} = {latex(total)}")
        C.append(row)
    w.m("C = A·B", C, caption="Результат", ans=True)
    return C


def transpose_op(A, w):
    w.h("Транспонирование матрицы: C = Aᵀ", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}")
    B = transpose(A)
    w.p("При транспонировании строки становятся столбцами: каждый элемент c_ij = a_ji.")
    w.m("C = Aᵀ", B, caption=f"C — матрица {n}×{m}", ans=True)
    return B


def rank(A, w):
    w.h("Ранг матрицы (метод Гаусса)", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}")
    w.p("Приводим матрицу к ступенчатому виду элементарными преобразованиями строк. "
        "Ранг равен числу ненулевых строк в ступенчатом виде.")
    M = [r[:] for r in A]
    cur = 0
    for c in range(n):
        piv = None
        for r in range(cur, m):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            continue
        if piv != cur:
            M[piv], M[cur] = M[cur], M[piv]
            w.p(f"Меняем местами строки {piv + 1} и {cur + 1}.")
            w.m("Матрица после перестановки", M)
        w.p(f"Ведущий элемент в столбце {c + 1}: a{cur + 1}{c + 1} = {latex(M[cur][c])}.")
        for r in range(cur + 1, m):
            if M[r][c] == 0:
                continue
            mc = M[r][c] / M[cur][c]
            w.l(f"R_{{{r + 1}}} \\leftarrow R_{{{r + 1}}} - \\left({latex(mc)}\\right) R_{{{cur + 1}}}")
            M[r] = [M[r][j] - mc * M[cur][j] for j in range(n)]
            w.m("Матрица после преобразования", M)
        cur += 1
    w.p(f"В ступенчатом виде {cur} ненулевых строк  →  ранг A = {cur}", ans=True)
    return cur