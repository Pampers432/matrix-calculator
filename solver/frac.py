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


def cell(f):
    n, d = f.numerator, f.denominator
    if d == 1:
        return f'<span class="val">{n}</span>'
    return f'<span class="frac"><span class="num">{n}</span><span class="den">{d}</span></span>'


def join_terms(vals):
    out = ""
    for i, v in enumerate(vals):
        if v >= 0:
            out += ("" if i == 0 else " + ") + numstr(v)
        else:
            out += ("" if i == 0 else " + ") + "(" + numstr(v) + ")"
    return out if out else "0"