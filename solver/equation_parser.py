from fractions import Fraction
import re

from .frac import parse_num


def is_identifier(value):
    if not isinstance(value, str) or not value:
        return False
    if not value[0].isalpha():
        return False
    return all(ch.isalnum() or ch == "_" for ch in value[1:])


def _tokenize(text):
    source = str(text).replace("−", "-").replace("·", "*").replace("×", "*").replace("•", "*")
    tokens = []
    i = 0
    while i < len(source):
        ch = source[i]
        if ch.isspace():
            i += 1
            continue
        if ch in "=+-*":
            tokens.append((ch, ch, i))
            i += 1
            continue
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
            tokens.append(("number", source[start:i], start))
            continue
        if ch.isalpha():
            start = i
            i += 1
            while i < len(source) and (source[i].isalnum() or source[i] == "_"):
                i += 1
            tokens.append(("name", source[start:i], start))
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
    _split_top_level(text[1:-1], ",")
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


def parse_matrix_equation(value):
    text = str(value).replace("−", "-").replace("·", "*").replace("×", "*")
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
    unknown = unknown.strip()
    if not is_identifier(unknown):
        raise ValueError("Неизвестная матрица должна иметь имя, начинающееся с буквы.")
    A = _parse_matrix_literal(matrix_text)
    B = _parse_matrix_literal(sides[1])
    if coefficient != 1:
        A = [[coefficient * value for value in row] for row in A]
    return {"A": A, "B": B, "unknown": unknown}


def _parse_equation_text(text):
    raw = str(text).strip()
    if not raw:
        raise ValueError("Пустое уравнение.")
    tokens = _tokenize(raw)
    equals = [i for i, token in enumerate(tokens) if token[0] == "="]
    if len(equals) != 1:
        raise ValueError("Уравнение должно содержать ровно один знак =.")
    split = equals[0]
    lhs, rhs = tokens[:split], tokens[split + 1:]
    if len(rhs) != 1 or rhs[0][0] != "name":
        position = rhs[0][2] if rhs else len(raw)
        raise ValueError(f"Справа от = должно быть имя матрицы (позиция {position + 1}).")

    terms = []
    i = 0
    sign = 1
    while i < len(lhs):
        if lhs[i][0] in ("+", "-"):
            sign = 1 if lhs[i][0] == "+" else -1
            i += 1
            if i >= len(lhs):
                raise ValueError("После знака ожидается член уравнения.")
        coefficient = Fraction(sign)
        if i < len(lhs) and lhs[i][0] == "number":
            number, i = _read_number(lhs, i, raw)
            coefficient *= number
            if i < len(lhs) and lhs[i][0] == "*":
                i += 1
        if i >= len(lhs):
            raise ValueError("После числового коэффициента ожидается имя матрицы.")
        first, i = _read_name(lhs, i, raw)
        second = None
        if i < len(lhs) and lhs[i][0] == "*":
            i += 1
            second, i = _read_name(lhs, i, raw)
        terms.append({"coefficient": coefficient, "first": first, "second": second})
        if i < len(lhs):
            if lhs[i][0] not in ("+", "-"):
                raise ValueError(
                    f"Между членами уравнения ожидался + или - (позиция {lhs[i][2] + 1})."
                )
            sign = 1
    if not terms and lhs:
        raise ValueError("Левая часть уравнения должна содержать матричный член.")
    return {"terms": terms, "rhs": rhs[0][1]}


def _normalise_structured(item):
    if not isinstance(item, dict):
        raise ValueError("Каждое уравнение должно быть строкой или объектом.")
    if item.get("text") is not None:
        return _parse_equation_text(item["text"])
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
        out.append({"coefficient": coefficient * sign, "first": first, "second": second})
    return {"terms": out, "rhs": rhs}


def parse_equations(value):
    if value is None:
        raise ValueError("Список уравнений не задан.")
    if isinstance(value, str):
        parts = [part.strip() for part in value.replace(";", "\n").splitlines()]
        return [_parse_equation_text(part) for part in parts if part]
    if not isinstance(value, list) or not value:
        raise ValueError("Нужна хотя бы одна система или одно уравнение.")
    result = []
    for item in value:
        if isinstance(item, str):
            parts = [part.strip() for part in item.replace(";", "\n").splitlines()]
            result.extend(_parse_equation_text(part) for part in parts if part)
        elif isinstance(item, dict):
            if item.get("text") is not None:
                parts = [part.strip() for part in str(item["text"]).replace(";", "\n").splitlines()]
                result.extend(_parse_equation_text(part) for part in parts if part)
            else:
                result.append(_normalise_structured(item))
        else:
            raise ValueError("Каждый элемент системы должен быть уравнением.")
    if not result:
        raise ValueError("Нужна хотя бы одна система или одно уравнение.")
    return result
