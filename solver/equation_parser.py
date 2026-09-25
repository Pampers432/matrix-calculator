from fractions import Fraction
import re

from .frac import latex, parse_num

MAX_DEPTH = 8

_REPLACEMENTS = (
    ("\u2212", "-"),
    ("\u00b7", "*"),
    ("\u00d7", "*"),
    ("\u2022", "*"),
    ("\u2019", "'"),
    ("\u02b9", "'"),
)
_SUPERSCRIPT_T = ("\u1d40", "\u1d57")
_TRANSVERSE_SUFFIXES = ("^T", "^t", "'") + _SUPERSCRIPT_T


def is_identifier(value):
    if not isinstance(value, str) or not value:
        return False
    if not value[0].isalpha():
        return False
    return all(ch.isalnum() or ch == "_" for ch in value[1:])


def name_latex(name):
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", name):
        return name
    return r"\mathrm{" + name.replace("_", r"\_") + "}"


def factor_text(factor):
    if factor["kind"] == "name":
        return name_latex(factor["name"]) + (r"^{T}" if factor["transposed"] else "")
    body = r"\left(" + _expression_latex(factor["terms"]) + r"\right)"
    return body + (r"^{T}" if factor["transposed"] else "")


def _expression_latex(terms):
    parts = []
    for term in terms:
        coefficient = term["coefficient"]
        body = term["body"]
        if abs(coefficient) != 1:
            body = latex(abs(coefficient)) + r"\," + body
        if not parts:
            parts.append(("-" if coefficient < 0 else "") + body)
        else:
            parts.append((" + " if coefficient > 0 else " - ") + body)
    return "".join(parts) if parts else "0"


def _annotate(terms):
    for term in terms:
        term["body"] = r"\cdot ".join(factor_text(f) for f in term["factors"])
    return terms


def _tokenize(text):
    source = str(text)
    for old, new in _REPLACEMENTS:
        source = source.replace(old, new)
    tokens = []
    i = 0
    while i < len(source):
        ch = source[i]
        if ch.isspace():
            i += 1
            continue
        if ch in "=+-*()":
            tokens.append((ch, ch, i, False))
            i += 1
            continue
        if ch == "'" or ch in _SUPERSCRIPT_T:
            tokens.append(("transpose", "T", i, False))
            i += 1
            continue
        if ch == "^":
            if i + 1 < len(source) and source[i + 1] in "Tt":
                tokens.append(("transpose", "T", i, False))
                i += 2
                continue
            raise ValueError(f"После ^ ожидался признак транспонирования T (позиция {i + 1}).")
        if ch.isdigit() or ch == ".":
            start = i
            while i < len(source) and (source[i].isdigit() or source[i] == "."):
                i += 1
            if i < len(source) and source[i] == "/":
                i += 1
                slash = i
                while i < len(source) and (source[i].isdigit() or source[i] == "."):
                    i += 1
                if i == slash:
                    raise ValueError(f"После / ожидалось число (позиция {i + 1}).")
            tokens.append(("number", source[start:i], start, False))
            continue
        if ch.isalpha():
            start = i
            i += 1
            while i < len(source) and source[i] not in _SUPERSCRIPT_T and (
                source[i].isalnum() or source[i] == "_"
            ):
                i += 1
            word = source[start:i]
            if len(word) == 1 and word in "Tt":
                tokens.append(("name", word, start, "bare"))
            else:
                tokens.append(("name", word, start, len(word) > 1 and word[-1] in "Tt"))
            continue
        raise ValueError(f"Недопустимый символ «{ch}» (позиция {i + 1}).")
    return tokens


def _read_name(tokens, index, text):
    if index >= len(tokens) or tokens[index][0] != "name":
        position = tokens[index][2] if index < len(tokens) else len(text)
        raise ValueError(f"Ожидалось имя матрицы (позиция {position + 1}).")
    return tokens[index][1], index + 1


def _read_number(tokens, index, text):
    if index >= len(tokens) or tokens[index][0] != "number":
        position = tokens[index][2] if index < len(tokens) else len(text)
        raise ValueError(f"Ожидалось числовое значение (позиция {position + 1}).")
    try:
        return parse_num(tokens[index][1]), index + 1
    except Exception as exc:
        raise ValueError(f"Неверное числовое значение: {exc}") from exc


def _read_postfix(tokens, index):
    if index < len(tokens):
        kind, _value, _pos, flag = tokens[index]
        if kind == "transpose" or (kind == "name" and flag == "bare"):
            return True, index + 1
    return False, index


def _flipped(factor):
    copy = dict(factor)
    copy["transposed"] = not factor["transposed"]
    return copy


def _parse_factor(tokens, index, text, names, depth):
    if index < len(tokens) and tokens[index][0] == "(":
        if depth >= MAX_DEPTH:
            raise ValueError("Слишком глубокая вложенность скобок.")
        index += 1
        inner, index = _parse_expression(tokens, index, text, names, depth + 1)
        if index >= len(tokens) or tokens[index][0] != ")":
            raise ValueError("В выражении не закрыта скобка.")
        index += 1
        transposed, index = _read_postfix(tokens, index)
        if len(inner) == 1:
            factors = [dict(factor) for factor in inner[0]["factors"]]
            if transposed:
                factors = [_flipped(factor) for factor in reversed(factors)]
            return factors, index
        return [{"kind": "expr", "terms": _annotate(inner), "transposed": transposed}], index
    if index < len(tokens) and tokens[index][0] == "number":
        raise ValueError(f"Число должно умножаться на матрицу (позиция {tokens[index][2] + 1}).")
    name, index = _read_name(tokens, index, text)
    transposed = False
    if len(name) > 1 and name[-1] in "Tt" and (names is None or name not in names):
        name = name[:-1]
        transposed = True
    extra, index = _read_postfix(tokens, index)
    return [{"kind": "name", "name": name, "transposed": transposed or extra}], index


def _parse_term(tokens, index, text, names, depth):
    coefficient = Fraction(1)
    if index < len(tokens) and tokens[index][0] == "number":
        coefficient, index = _read_number(tokens, index, text)
        if index < len(tokens) and tokens[index][0] == "*":
            index += 1
    factors, index = _parse_factor(tokens, index, text, names, depth)
    while index < len(tokens) and tokens[index][0] == "*":
        index += 1
        more, index = _parse_factor(tokens, index, text, names, depth)
        factors.extend(more)
    return {"coefficient": coefficient, "factors": factors}, index


def _parse_expression(tokens, index, text, names, depth=0):
    terms = []
    sign = 1
    while True:
        term, index = _parse_term(tokens, index, text, names, depth)
        term["coefficient"] = sign * term["coefficient"]
        terms.append(term)
        if index < len(tokens) and tokens[index][0] in ("+", "-"):
            sign = 1 if tokens[index][0] == "+" else -1
            index += 1
            continue
        break
    return terms, index


def _split_top_level(text, separators):
    stack = []
    parts = []
    start = 0
    closing = {"]": "[", ")": "("}
    for index, char in enumerate(text):
        if char in "[(":
            stack.append(char)
        elif char in "])":
            if not stack or stack[-1] != closing[char]:
                raise ValueError(f"Неверная скобка в позиции {index + 1}.")
            stack.pop()
        elif char in separators and not stack:
            parts.append(text[start:index].strip())
            start = index + 1
    if stack:
        raise ValueError("В выражении не закрыта скобка.")
    parts.append(text[start:].strip())
    return parts


def _is_wrapped(text):
    if len(text) < 2:
        return False
    return (text[0] == "[" and text[-1] == "]") or (text[0] == "(" and text[-1] == ")")


def _parse_matrix_row(text):
    text = text.strip()
    if _is_wrapped(text):
        text = text[1:-1].strip()
    if not text:
        raise ValueError("В матрице не может быть пустой строки.")
    values = []
    for part in re.split(r"\s*,\s*", text):
        value = re.sub(r"\s+", "", part)
        if not value:
            raise ValueError("В матрице не может быть пустого элемента.")
        try:
            values.append(parse_num(value))
        except Exception as exc:
            raise ValueError(f"Неверный элемент матрицы: {exc}") from exc
    return values


def _parse_matrix_literal(text):
    text = text.strip()
    if not _is_wrapped(text):
        raise ValueError("Матрица должна быть записана в квадратных скобках: [[1,2],[3,4]].")
    inner = text[1:-1].strip()
    if not inner:
        raise ValueError("Матрица не может быть пустой.")
    first = inner[0]
    separators = ",;\n" if first in "[(" else ";\n"
    rows = []
    for row_text in _split_top_level(inner, separators):
        rows.append(_parse_matrix_row(row_text))
    if not rows or any(not row for row in rows):
        raise ValueError("Матрица не может быть пустой.")
    if len({len(row) for row in rows}) != 1:
        raise ValueError("Все строки матрицы должны быть одинаковой длины.")
    return rows


def _split_transpose(text):
    stripped = text.strip()
    for suffix in _TRANSVERSE_SUFFIXES:
        if len(stripped) > len(suffix) and stripped.endswith(suffix):
            return stripped[: -len(suffix)].strip(), True
    if len(stripped) > 1 and stripped[-1] in "Tt" and stripped[-2] in "])":
        return stripped[:-1].strip(), True
    return stripped, False


def parse_matrix_equation(value):
    text = str(value)
    for old, new in _REPLACEMENTS:
        text = text.replace(old, new)
    text = text.strip()
    if not text:
        raise ValueError("Введите матричное уравнение.")
    sides = _split_top_level(text, "=")
    if len(sides) != 2 or not sides[0] or not sides[1]:
        raise ValueError("Уравнение должно содержать ровно один знак =.")
    parts = _split_top_level(sides[0], "*")
    if len(parts) == 2:
        coefficient = Fraction(1)
        matrix_text, unknown = parts
    elif len(parts) == 3:
        try:
            coefficient = parse_num(parts[0])
        except Exception as exc:
            raise ValueError(f"Неверный коэффициент матрицы: {exc}") from exc
        matrix_text, unknown = parts[1], parts[2]
    else:
        raise ValueError("Левая часть должна иметь вид A*X или k*A*X.")
    unknown, unknown_t = _split_transpose(unknown)
    if not is_identifier(unknown):
        raise ValueError("Неизвестная матрица должна иметь имя, начинающееся с буквы.")
    matrix_text, matrix_t = _split_transpose(matrix_text)
    rhs_text, rhs_t = _split_transpose(sides[1])
    A = _parse_matrix_literal(matrix_text)
    B = _parse_matrix_literal(rhs_text)
    if matrix_t:
        A = [list(row) for row in zip(*A)]
    if rhs_t:
        B = [list(row) for row in zip(*B)]
    if coefficient != 1:
        A = [[coefficient * value for value in row] for row in A]
    return {"A": A, "B": B, "unknown": unknown, "unknown_transposed": unknown_t}


def _parse_equation_text(text, names=None):
    raw = str(text).strip()
    if not raw:
        raise ValueError("Пустое уравнение.")
    tokens = _tokenize(raw)
    equals = [i for i, token in enumerate(tokens) if token[0] == "="]
    if len(equals) != 1:
        raise ValueError("Уравнение должно содержать ровно один знак =.")
    split = equals[0]
    lhs_tokens, rhs_tokens = tokens[:split], tokens[split + 1:]
    if not lhs_tokens:
        raise ValueError("Левая часть уравнения должна содержать матричный член.")
    if not rhs_tokens:
        raise ValueError("Справа от = должно быть выражение из известных матриц.")
    lhs, index = _parse_expression(lhs_tokens, 0, raw, names)
    if index != len(lhs_tokens):
        raise ValueError(f"Лишние символы в левой части (позиция {lhs_tokens[index][2] + 1}).")
    rhs, index = _parse_expression(rhs_tokens, 0, raw, names)
    if index != len(rhs_tokens):
        raise ValueError(f"Лишние символы в правой части (позиция {rhs_tokens[index][2] + 1}).")
    equation = {"terms": _annotate(lhs), "rhs_terms": _annotate(rhs)}
    _describe(equation)
    return equation


def _describe(equation):
    for group in (equation["terms"], equation["rhs_terms"]):
        for term in group:
            refs = [f for f in term["factors"] if f["kind"] == "name"]
            term["first"] = refs[0]["name"] if refs else None
            term["second"] = refs[1]["name"] if len(refs) > 1 else None
            term["first_transposed"] = bool(refs and refs[0]["transposed"])
            term["second_transposed"] = bool(len(refs) > 1 and refs[1]["transposed"])
    equation["rhs"] = None
    equation["rhs_transposed"] = False
    single = equation["rhs_terms"]
    if len(single) == 1 and single[0]["coefficient"] == 1 and len(single[0]["factors"]) == 1:
        factor = single[0]["factors"][0]
        if factor["kind"] == "name":
            equation["rhs"] = factor["name"]
            equation["rhs_transposed"] = factor["transposed"]
    equation["rhs_latex"] = _expression_latex(equation["rhs_terms"])
    equation["latex"] = _expression_latex(equation["terms"]) + " = " + equation["rhs_latex"]


def _normalise_structured(item, names=None):
    if not isinstance(item, dict):
        raise ValueError("Каждое уравнение должно быть строкой или объектом.")
    if item.get("text") is not None:
        return _parse_equation_text(item["text"], names)
    terms = item.get("terms")
    rhs = item.get("rhs")
    if not isinstance(terms, list) or not isinstance(rhs, str) or not is_identifier(rhs):
        raise ValueError("Структурное уравнение должно содержать terms и имя матрицы rhs.")
    out = []
    for term in terms:
        if not isinstance(term, dict):
            raise ValueError("Каждый член уравнения должен быть объектом.")
        first = term.get("first", term.get("left"))
        second = term.get("second", term.get("right"))
        known_name = term.get("known")
        unknown_name = term.get("unknown")
        if first is None and known_name is not None:
            first = known_name
            if second is None:
                second = unknown_name
        elif first is None and unknown_name is not None:
            first = unknown_name
        if second is None and unknown_name is not None and first != unknown_name:
            second = unknown_name
        if not is_identifier(first) or (second is not None and not is_identifier(second)):
            raise ValueError("В члене уравнения указано неверное имя матрицы.")
        try:
            coefficient = parse_num(term.get("coefficient", 1))
        except Exception as exc:
            raise ValueError(f"Неверный коэффициент члена: {exc}") from exc
        sign = term.get("sign", term.get("op", 1))
        if sign in ("+", 1, "1", "+1"):
            sign = 1
        elif sign in ("-", -1, "-1"):
            sign = -1
        else:
            raise ValueError("Знак члена должен быть + или -.")
        factors = [{"kind": "name", "name": first, "transposed": bool(term.get("transposed"))}]
        if second is not None:
            factors.append({
                "kind": "name",
                "name": second,
                "transposed": bool(term.get("second_transposed")),
            })
        out.append({"coefficient": coefficient * sign, "factors": factors})
    rhs_terms = [{
        "coefficient": Fraction(1),
        "factors": [{"kind": "name", "name": rhs, "transposed": bool(item.get("rhs_transposed"))}],
    }]
    equation = {"terms": _annotate(out), "rhs_terms": _annotate(rhs_terms)}
    _describe(equation)
    return equation


def parse_equations(value, names=None):
    if value is None:
        raise ValueError("Список уравнений не задан.")
    if isinstance(value, str):
        parts = [part.strip() for part in value.replace(";", "\n").splitlines()]
        return [_parse_equation_text(part, names) for part in parts if part]
    if not isinstance(value, list) or not value:
        raise ValueError("Нужна хотя бы одна система или одно уравнение.")
    result = []
    for item in value:
        if isinstance(item, str):
            parts = [part.strip() for part in item.replace(";", "\n").splitlines()]
            result.extend(_parse_equation_text(part, names) for part in parts if part)
        elif isinstance(item, dict):
            if item.get("text") is not None:
                parts = [part.strip() for part in str(item["text"]).replace(";", "\n").splitlines()]
                result.extend(_parse_equation_text(part, names) for part in parts if part)
            else:
                result.append(_normalise_structured(item, names))
        else:
            raise ValueError("Каждый элемент системы должен быть уравнением.")
    if not result:
        raise ValueError("Нужна хотя бы одна система или одно уравнение.")
    return result
