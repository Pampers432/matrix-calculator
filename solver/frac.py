from fractions import Fraction


def parse_num(v):
    if isinstance(v, Fraction):
        return v
    if isinstance(v, bool):
        return Fraction(int(v))
    if isinstance(v, (int, float)):
        return Fraction(str(v))
    s = str(v).strip()
    if not s:
        return Fraction(0)
    if "/" in s:
        a, b = s.split("/", 1)
        den = parse_num(b)
        if den == 0:
            raise ValueError("Деление на ноль в значении матрицы")
        return parse_num(a) / den
    return Fraction(s)


def numstr(f):
    if f is None:
        return ""
    if f.denominator == 1:
        return str(f.numerator)
    if f.numerator == 0:
        return "0"
    return f"{f.numerator}/{f.denominator}"


def latex(f):
    n, d = f.numerator, f.denominator
    if d == 1:
        return str(n)
    return f"\\frac{{{n}}}{{{d}}}"


def ljoin(vals):
    parts = []
    for i, v in enumerate(vals):
        if v >= 0:
            parts.append((" + " if i else "") + latex(v))
        else:
            parts.append((" - " if i else "-") + latex(-v))
    return "".join(parts) if parts else "0"