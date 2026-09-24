from fractions import Fraction

from .frac import numstr, join_terms
from .matutil import transpose


def _opfmt(v):
    return "(" + numstr(v) + ")" if v < 0 else numstr(v)


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
            sign = " + " if op == "+" else " − "
            w.p(f"c{i + 1}{j + 1} = a{i + 1}{j + 1}{sign}b{i + 1}{j + 1} = "
                f"{_opfmt(A[i][j])}{sign}{_opfmt(B[i][j])} = {numstr(v)}")
        C.append(row)
    w.m(f"C = A {op} B", C, caption="Результат", ans=True)
    return C


def scalar_mul(k, A, w):
    w.h(f"Умножение матрицы на число: C = {numstr(k)}·A", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}")
    w.p("Каждый элемент матрицы умножается на число k:")
    C = []
    for i in range(m):
        row = []
        for j in range(n):
            v = k * A[i][j]
            row.append(v)
            w.p(f"c{i + 1}{j + 1} = {numstr(k)}·a{i + 1}{j + 1} = {numstr(k)}·{_opfmt(A[i][j])} = {numstr(v)}")
        C.append(row)
    w.m("C = k·A", C, caption="Результат", ans=True)
    return C


def mat_mul(A, B, w):
    w.h("Умножение матриц: C = A·B", 1)
    m, p = len(A), len(A[0])
    n = len(B[0])
    w.m("A", A, caption=f"A — матрица {m}×{p}")
    w.m("B", B, caption=f"B — матрица {p}×{n}")
    w.p(f"Число столбцов A ({p}) равно числу строк B ({p}) — матрицы согласованы, "
        f"результат C имеет размер {m}×{n}.")
    w.p("Элемент c_ij = сумма произведений элементов i-й строки A на j-й столбец B.")
    C = []
    for i in range(m):
        row = []
        for j in range(n):
            terms = [A[i][k] * B[k][j] for k in range(p)]
            total = sum(terms, Fraction(0))
            row.append(total)
            formulas = " + ".join(f"a{i + 1}{k + 1}·b{k + 1}{j + 1}" for k in range(p))
            w.p(f"c{i + 1}{j + 1} = {formulas} = {join_terms(terms)} = {numstr(total)}")
        C.append(row)
    w.m("C = A·B", C, caption="Результат", ans=True)
    return C


def transpose_op(A, w):
    w.h("Транспонирование матрицы: C = Aᵀ", 1)
    m, n = len(A), len(A[0])
    w.m("A", A, caption=f"A — матрица {m}×{n}")
    B = transpose(A)
    w.p("При транспонировании строки становятся столбцами: c_ij = a_ji.")
    if m * n <= 16:
        for i in range(n):
            for j in range(m):
                w.p(f"c{i + 1}{j + 1} = a{j + 1}{i + 1} = {numstr(A[j][i])}")
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
        w.p(f"Ведущий элемент в столбце {c + 1}: a{cur + 1}{c + 1} = {numstr(M[cur][c])}")
        for r in range(cur + 1, m):
            if M[r][c] == 0:
                continue
            mc = M[r][c] / M[cur][c]
            w.p(f"R{r + 1} ← R{r + 1} − ({numstr(mc)})·R{cur + 1}")
            M[r] = [M[r][j] - mc * M[cur][j] for j in range(n)]
            w.m("Матрица после преобразования", M)
        cur += 1
    w.p(f"В ступенчатом виде {cur} ненулевых строк  →  ранг A = {cur}", ans=True)
    return cur