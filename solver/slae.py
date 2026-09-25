from fractions import Fraction

from .frac import latex
from .det import laplace


def _a_row(aug):
    n = len(aug[0]) - 1
    return [aug[i][:n] + [{"sep": True}] + aug[i][n:] for i in range(len(aug))]


def _b_col(b):
    return [[x] for x in b]


def cramer(A, b, w):
    n = len(A)
    w.h("Решение системы линейных уравнений (СЛАУ) методом Крамера", 1)
    w.m("A", A, caption="Матрица системы (n×n)")
    w.m("b", _b_col(b), caption="Столбец свободных членов")

    w.h("Главный определитель системы Δ = det(A)", 2)
    D = laplace(A, w)
    w.l(f"\\Delta = {latex(D)}", ans=True)
    if D == 0:
        w.p("Δ = 0 — либо система несовместна, либо имеет бесконечно много решений. "
            "Метод Крамера неприменим, воспользуйтесь методом Гаусса.", ans=True)
        return

    xs = []
    for i in range(n):
        w.h(f"Определитель Δ{i + 1}: заменяем {i + 1}-й столбец A на столбец b", 2)
        Di = [r[:] for r in A]
        for r in range(n):
            Di[r][i] = b[r]
        w.m(f"\\Delta_{{{i + 1}}}", Di)
        d = laplace(Di, w, label=f"\\Delta_{{{i + 1}}}", heading=False)
        w.l(f"\\Delta_{{{i + 1}}} = {latex(d)}", ans=True)
        x = d / D
        xs.append(x)
        w.l(f"x_{{{i + 1}}} = \\dfrac{{\\Delta_{{{i + 1}}}}}{{\\Delta}} = \\dfrac{{{latex(d)}}}{{{latex(D)}}} = {latex(x)}", ans=True)

    w.h("Решение системы", 1)
    w.m("X", _b_col(xs), caption="Вектор неизвестных", ans=True)
    return xs


def gauss(A, b, w):
    m = len(A)
    nc = len(A[0])
    w.h("Решение системы линейных уравнений (СЛАУ) методом Гаусса", 1)
    w.m("A", A, caption=f"Матрица системы ({m}×{nc}, {m} уравнений, {nc} неизвестных)")
    w.m("b", _b_col(b), caption="Столбец свободных членов")
    aug = [A[i][:] + [b[i]] for i in range(m)]
    w.p("Составляем расширенную матрицу [A | b]:")
    w.m("[A | b]", _a_row(aug))

    pivots = []
    cur = 0
    for c in range(nc):
        piv = None
        for r in range(cur, m):
            if aug[r][c] != 0:
                piv = r
                break
        if piv is None:
            w.p(f"В столбце {c + 1} среди оставшихся строк ненулевых нет → x{c + 1} будет свободной переменной.")
            continue
        if piv != cur:
            aug[piv], aug[cur] = aug[cur], aug[piv]
            w.p(f"Меняем местами строки {piv + 1} и {cur + 1}.")
            w.m("[A | b]", _a_row(aug))
        w.p(f"Нормируем строку {cur + 1}: делим на ведущий элемент {latex(aug[cur][c])}.")
        aug[cur] = [v / aug[cur][c] for v in aug[cur]]
        w.m("[A | b]", _a_row(aug))
        for r in range(m):
            if r == cur or aug[r][c] == 0:
                continue
            mm = aug[r][c]
            w.l(f"R_{{{r + 1}}} \\leftarrow R_{{{r + 1}}} - \\left({latex(mm)}\\right) R_{{{cur + 1}}}")
            aug[r] = [aug[r][j] - mm * aug[cur][j] for j in range(nc + 1)]
            w.m("[A | b]", _a_row(aug))
        pivots.append(c)
        cur += 1

    for r in range(cur, m):
        if all(aug[r][j] == 0 for j in range(nc)) and aug[r][nc] != 0:
            w.l(f"0 = {latex(aug[r][nc])} \\quad \\text{{— противоречие}}")
            w.p("Система несовместна — решений нет.", ans=True)
            return

    r = len(pivots)
    w.p(f"Ранг матрицы системы равен рангу расширенной матрицы и равен {r}.")

    if r == nc:
        w.p("Все неизвестные главные — система имеет единственное решение:")
        X = [None] * nc
        for i, c in enumerate(pivots):
            X[c] = aug[i][nc]
        for i in range(nc):
            w.l(f"x_{{{i + 1}}} = {latex(X[i])}")
        w.m("X", _b_col(X), caption="Вектор неизвестных", ans=True)
        return

    free = [c for c in range(nc) if c not in pivots]
    params = ["t_1", "t_2", "t_3", "t_4", "t_5", "t_6"]
    names = {c: params[idx] for idx, c in enumerate(free)}
    w.p("Есть свободные переменные — система имеет бесконечно много решений. "
        "Полагаем " + ", ".join(f"x{c + 1} = t_{{{idx + 1}}}" for idx, c in enumerate(free)) + ".")
    sol = {}
    for i in reversed(range(len(pivots))):
        c = pivots[i]
        parts = []
        const = aug[i][nc]
        if const != 0:
            parts.append(latex(const))
        for j in range(nc):
            if j == c or aug[i][j] == 0:
                continue
            t = -aug[i][j]
            nm = names[j]
            cs = nm if abs(t) == 1 else f"{latex(abs(t))}\\,{nm}"
            if not parts:
                parts.append(cs if t > 0 else "-" + cs)
            else:
                parts.append((" + " if t > 0 else " - ") + cs)
        sol[c] = "".join(parts) if parts else "0"
        w.l(f"x_{{{c + 1}}} = {sol[c]}")

    w.h("Решение (в параметрическом виде)", 1)
    rows = []
    for i in range(nc):
        if i in sol:
            rows.append([sol[i]])
        else:
            rows.append([names[i]])
    w.mraw("X", rows, caption="Вектор неизвестных (t — произвольные параметры)", ans=True)