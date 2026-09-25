import re
from fractions import Fraction

from .equation_parser import is_identifier, parse_equations
from .frac import latex, parse_num
from .matutil import transpose


MAX_MATRIX_SIZE = 6
MAX_SCALAR_EQUATIONS = 144
MAX_SCALAR_VARIABLES = 144


def _shape(M):
    return len(M), len(M[0])


def _parse_matrix(data):
    if not isinstance(data, list) or not data:
        raise ValueError("Матрица не задана.")
    matrix = []
    for row in data:
        if not isinstance(row, list):
            row = [row]
        if not row:
            raise ValueError("Матрица не может содержать пустые строки.")
        matrix.append([parse_num(value) for value in row])
    if len(matrix) > MAX_MATRIX_SIZE or any(len(row) > MAX_MATRIX_SIZE for row in matrix):
        raise ValueError(f"Размер матриц ограничен {MAX_MATRIX_SIZE}×{MAX_MATRIX_SIZE}.")
    if len({len(row) for row in matrix}) > 1:
        raise ValueError("Все строки матрицы должны быть одинаковой длины.")
    return matrix


def _positive_int(value, label):
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"{label} должно быть положительным целым числом.")
    try:
        number = int(value)
    except Exception as exc:
        raise ValueError(f"{label} должно быть положительным целым числом.") from exc
    if number < 1 or number > MAX_MATRIX_SIZE or str(value).strip() != str(number):
        raise ValueError(f"{label} должно быть от 1 до {MAX_MATRIX_SIZE}.")
    return number


def _parse_unknowns(value):
    if value is None:
        return {}
    result = {}
    if isinstance(value, dict):
        items = []
        for name, spec in value.items():
            if isinstance(spec, dict):
                item = dict(spec)
                item.setdefault("name", name)
            elif spec is None:
                item = {"name": name}
            else:
                item = {"name": name, "rows": spec, "cols": spec}
            items.append(item)
    elif isinstance(value, list):
        items = value
    else:
        raise ValueError("Неизвестные матрицы должны быть заданы списком или объектом.")
    for item in items:
        if isinstance(item, str):
            name, rows, cols = item, None, None
        elif isinstance(item, dict):
            name = item.get("name")
            size = item.get("size")
            if isinstance(size, (list, tuple)) and len(size) == 2:
                rows, cols = size
            else:
                rows, cols = item.get("rows"), item.get("cols")
        else:
            raise ValueError("Описание неизвестной матрицы имеет неверный формат.")
        if not is_identifier(name):
            raise ValueError("Имя неизвестной матрицы должно начинаться с буквы и содержать только буквы, цифры и _.")
        if name in result:
            raise ValueError(f"Неизвестная матрица {name} объявлена несколько раз.")
        result[name] = (
            _positive_int(rows, f"Число строк {name}"),
            _positive_int(cols, f"Число столбцов {name}"),
        )
    return result


def _parse_known(value):
    if value is None:
        raise ValueError("Известные матрицы не заданы.")
    if isinstance(value, dict):
        entries = list(value.items())
    elif isinstance(value, list) and value:
        entries = []
        for item in value:
            if not isinstance(item, dict) or "name" not in item:
                raise ValueError("Каждая известная матрица должна иметь имя и значение data.")
            if "data" not in item and "matrix" not in item:
                raise ValueError(f"Для известной матрицы {item.get('name')} не задано значение data.")
            entries.append((item["name"], item.get("data", item.get("matrix"))))
    else:
        raise ValueError("Известные матрицы должны быть заданы объектом или списком.")
    if not entries:
        raise ValueError("Известные матрицы должны быть заданы объектом или списком.")
    result = {}
    for name, data in entries:
        if not is_identifier(name):
            raise ValueError(f"Недопустимое имя матрицы: {name}.")
        if name in result:
            raise ValueError(f"Матрица {name} объявлена несколько раз.")
        result[name] = _parse_matrix(data)
    return result


def _infer_shape(shapes, name, rows, cols):
    current_rows, current_cols = shapes[name]
    if current_rows is not None and rows is not None and current_rows != rows:
        raise ValueError(f"Для матрицы {name} указаны разные числа строк.")
    if current_cols is not None and cols is not None and current_cols != cols:
        raise ValueError(f"Для матрицы {name} указаны разные числа столбцов.")
    shapes[name] = (
        current_rows if current_rows is not None else rows,
        current_cols if current_cols is not None else cols,
    )


def _prepare_equations(equations, known, shapes):
    collisions = set(known) & set(shapes)
    if collisions:
        raise ValueError(
            "Имя не может одновременно обозначать известную и неизвестную матрицу: "
            + ", ".join(sorted(collisions))
        )
    names = list(shapes)
    for equation in equations:
        if equation["rhs"] not in known:
            raise ValueError(f"Правая часть {equation['rhs']} должна быть известной матрицей.")
        for term in equation["terms"]:
            for name in (term["first"], term.get("second")):
                if name is not None and name not in known and name not in shapes:
                    shapes[name] = (None, None)
                    names.append(name)

    prepared = []
    for equation in equations:
        rhs = known[equation["rhs"]]
        rhs_rows, rhs_cols = _shape(rhs)
        terms = []
        for raw in equation["terms"]:
            first = raw["first"]
            second = raw.get("second")
            first_known = first in known
            second_known = second in known if second is not None else False
            first_unknown = first in shapes
            second_unknown = second in shapes if second is not None else False
            if first_known and second_known:
                raise ValueError("Произведение двух известных матриц не входит в линейную систему.")
            if first_unknown and second_unknown:
                raise ValueError("Произведение двух неизвестных матриц нелинейно и не поддерживается.")
            if first_known and second_unknown:
                kind = "known_unknown"
                matrix = known[first]
                unknown = second
                output_rows, output_cols = len(matrix), shapes[unknown][1] or rhs_cols
                _infer_shape(shapes, unknown, len(matrix[0]), rhs_cols)
            elif first_unknown and second_known:
                kind = "unknown_known"
                matrix = known[second]
                unknown = first
                output_rows, output_cols = shapes[unknown][0] or rhs_rows, len(matrix[0])
                _infer_shape(shapes, unknown, rhs_rows, len(matrix))
            elif first_unknown and second is None:
                kind = "unknown"
                unknown = first
                _infer_shape(shapes, unknown, rhs_rows, rhs_cols)
                output_rows, output_cols = shapes[unknown]
            elif first_known and second is None:
                kind = "known"
                output_rows, output_cols = _shape(known[first])
            else:
                raise ValueError("Член уравнения должен иметь вид A·X, X·A, A или X.")
            if (output_rows, output_cols) != (rhs_rows, rhs_cols):
                raise ValueError(
                    f"Размеры члена {first}"
                    + (f"·{second}" if second else "")
                    + f" ({output_rows}×{output_cols}) не совпадают с правой частью "
                    + f"({rhs_rows}×{rhs_cols})."
                )
            terms.append({
                "kind": kind,
                "coefficient": raw["coefficient"],
                "first": first,
                "second": second,
            })
        prepared.append({"rhs": rhs, "rhs_name": equation["rhs"], "terms": terms})

    for name, shape in shapes.items():
        if shape[0] is None or shape[1] is None:
            raise ValueError(f"Не удалось определить размеры неизвестной матрицы {name}.")
    return prepared, names


def _build_scalar_system(prepared, shapes, names, known):
    offsets = {}
    total = 0
    for name in names:
        offsets[name] = total
        total += shapes[name][0] * shapes[name][1]
    if total > MAX_SCALAR_VARIABLES:
        raise ValueError(f"Слишком много неизвестных элементов ({total}). Максимум: {MAX_SCALAR_VARIABLES}.")
    rows = []
    values = []
    for equation in prepared:
        rhs = equation["rhs"]
        rhs_rows, rhs_cols = _shape(rhs)
        for r in range(rhs_rows):
            for c in range(rhs_cols):
                row = [Fraction(0) for _ in range(total)]
                constant = Fraction(0)
                for term in equation["terms"]:
                    coefficient = term["coefficient"]
                    kind = term["kind"]
                    if kind == "known_unknown":
                        known_name = term["first"]
                        unknown_name = term["second"]
                        matrix = known[known_name]
                        unknown_shape = shapes[unknown_name]
                        offset = offsets[unknown_name]
                        for k in range(unknown_shape[0]):
                            index = offset + k * unknown_shape[1] + c
                            row[index] += coefficient * matrix[r][k]
                    elif kind == "unknown_known":
                        known_name = term["second"]
                        unknown_name = term["first"]
                        matrix = known[known_name]
                        unknown_shape = shapes[unknown_name]
                        offset = offsets[unknown_name]
                        for k in range(unknown_shape[1]):
                            index = offset + r * unknown_shape[1] + k
                            row[index] += coefficient * matrix[k][c]
                    elif kind == "unknown":
                        unknown_name = term["first"]
                        unknown_shape = shapes[unknown_name]
                        row[offsets[unknown_name] + r * unknown_shape[1] + c] += coefficient
                    else:
                        constant += coefficient * known[term["first"]][r][c]
                if len(rows) >= MAX_SCALAR_EQUATIONS:
                    raise ValueError(
                        f"Слишком много скалярных уравнений ({len(rows) + 1}). Максимум: {MAX_SCALAR_EQUATIONS}."
                    )
                rows.append(row)
                values.append(rhs[r][c] - constant)
    return rows, values


def _rref(C, D):
    m = len(C)
    n = len(C[0]) if m else 0
    q = len(D[0]) if m else 0
    matrix = [C[i][:] + D[i][:] for i in range(m)]
    pivots = []
    pivot_rows = {}
    current = 0
    for column in range(n):
        pivot = None
        for row in range(current, m):
            if matrix[row][column] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        if pivot != current:
            matrix[pivot], matrix[current] = matrix[current], matrix[pivot]
        divisor = matrix[current][column]
        matrix[current] = [value / divisor for value in matrix[current]]
        for row in range(m):
            if row == current or matrix[row][column] == 0:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                matrix[row][j] - multiplier * matrix[current][j]
                for j in range(n + q)
            ]
        pivots.append(column)
        pivot_rows[column] = current
        current += 1
    free = [column for column in range(n) if column not in pivot_rows]
    consistent = True
    for row in matrix:
        if all(value == 0 for value in row[:n]) and any(value != 0 for value in row[n:]):
            consistent = False
            break
    return matrix, pivots, pivot_rows, free, consistent


def _augmented_rows(C, D, separator):
    return [C[i][:] + [{"sep": True}] + D[i][:] for i in range(len(C))]


def _rref_with_steps(C, D, w, title):
    m = len(C)
    n = len(C[0]) if m else 0
    q = len(D[0]) if m else 0
    if m > 16 or n + q > 16:
        matrix, pivots, pivot_rows, free, consistent = _rref(C, D)
        w.m(title, [row[:n] + [{"sep": True}] + row[n:] for row in matrix])
        return matrix, pivots, pivot_rows, free, consistent
    w.m(title, _augmented_rows(C, D, True))
    matrix = [C[i][:] + D[i][:] for i in range(m)]
    pivots = []
    pivot_rows = {}
    current = 0
    for column in range(n):
        pivot = None
        for row in range(current, m):
            if matrix[row][column] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        if pivot != current:
            matrix[pivot], matrix[current] = matrix[current], matrix[pivot]
            w.p(f"Меняем местами строки {pivot + 1} и {current + 1}.")
            w.m(title, [row[:n] + [{"sep": True}] + row[n:] for row in matrix])
        divisor = matrix[current][column]
        w.p(
            f"Нормируем строку {current + 1} по ведущему элементу "
            f"{latex(matrix[current][column])}."
        )
        w.l(
            f"R_{{{current + 1}}} \\leftarrow "
            f"\\frac{{1}}{{{latex(matrix[current][column])}}}R_{{{current + 1}}}"
        )
        matrix[current] = [value / divisor for value in matrix[current]]
        w.m(title, [row[:n] + [{"sep": True}] + row[n:] for row in matrix])
        for row in range(m):
            if row == current or matrix[row][column] == 0:
                continue
            multiplier = matrix[row][column]
            w.l(
                f"R_{{{row + 1}}} \\leftarrow R_{{{row + 1}}} - "
                f"\\left({latex(multiplier)}\\right)R_{{{current + 1}}}"
            )
            matrix[row] = [
                matrix[row][j] - multiplier * matrix[current][j]
                for j in range(n + q)
            ]
            w.m(title, [row[:n] + [{"sep": True}] + row[n:] for row in matrix])
        pivots.append(column)
        pivot_rows[column] = current
        current += 1
    free = [column for column in range(n) if column not in pivot_rows]
    consistent = not any(
        all(value == 0 for value in row[:n]) and any(value != 0 for value in row[n:])
        for row in matrix
    )
    return matrix, pivots, pivot_rows, free, consistent


def _joined_augmented(matrix, n):
    return [row[:n] + [{"sep": True}] + row[n:] for row in matrix]


def _rref_by_addition(C, D, w, title):
    m = len(C)
    n = len(C[0]) if m else 0
    q = len(D[0]) if m else 0
    matrix = [C[i][:] + D[i][:] for i in range(m)]
    detailed = m <= 16 and n + q <= 16

    def show():
        if detailed:
            w.m(title, _joined_augmented(matrix, n))

    def add_equation(row, pivot, multiplier, column=None):
        matrix[row] = [
            matrix[row][j] + multiplier * matrix[pivot][j]
            for j in range(n + q)
        ]
        if not detailed:
            return
        if multiplier < 0:
            operation = (
                rf"\text{{уравнение }} {row + 1} \leftarrow "
                rf"\text{{уравнение }} {row + 1} - {latex(-multiplier)}\cdot"
                rf"\text{{уравнение }} {pivot + 1}"
            )
            action = (
                f"Вычитаем из уравнения {row + 1} уравнение {pivot + 1}, "
                f"умноженное на {latex(-multiplier)}"
            )
        else:
            operation = (
                rf"\text{{уравнение }} {row + 1} \leftarrow "
                rf"\text{{уравнение }} {row + 1} + {latex(multiplier)}\cdot"
                rf"\text{{уравнение }} {pivot + 1}"
            )
            action = (
                f"Прибавляем к уравнению {row + 1} уравнение {pivot + 1}, "
                f"умноженное на {latex(multiplier)}"
            )
        prefix = "Чтобы уравнять коэффициенты в столбце " + str(column + 1) + ", " if column is not None else ""
        w.p(prefix + action + ".")
        w.l(operation)
        show()

    show()
    pivots = []
    pivot_rows = {}
    current = 0
    for column in range(n):
        if current >= m:
            break
        if matrix[current][column] == 0:
            donor = next(
                (row for row in range(current + 1, m) if matrix[row][column] != 0),
                None,
            )
            if donor is None:
                continue
            add_equation(current, donor, Fraction(1))
        divisor = matrix[current][column]
        for row in range(current + 1, m):
            factor = matrix[row][column]
            if factor != 0:
                multiplier = -factor / divisor
                add_equation(row, current, multiplier, column)
        if divisor != 1:
            reciprocal = Fraction(1, 1) / divisor
            matrix[current] = [value * reciprocal for value in matrix[current]]
            if detailed:
                w.p(
                    f"Нормируем уравнение {current + 1} по ведущему элементу "
                    f"{latex(divisor)}."
                )
                w.l(
                    rf"\text{{уравнение }} {current + 1} \leftarrow "
                    rf"\frac{{1}}{{{latex(divisor)}}}\cdot"
                    rf"\text{{уравнение }} {current + 1}"
                )
                show()
        pivots.append(column)
        pivot_rows[column] = current
        current += 1

    for column in reversed(pivots):
        pivot_row = pivot_rows[column]
        for row in range(pivot_row):
            factor = matrix[row][column]
            if factor != 0:
                add_equation(row, pivot_row, -factor)

    if not detailed:
        w.m(title, _joined_augmented(matrix, n))
    free = [column for column in range(n) if column not in pivot_rows]
    consistent = not any(
        all(value == 0 for value in row[:n]) and any(value != 0 for value in row[n:])
        for row in matrix
    )
    return matrix, pivots, pivot_rows, free, consistent


def _name_latex(name):
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", name):
        return name
    escaped = name.replace("_", r"\_")
    return r"\mathrm{" + escaped + "}"


def _direct_term_key(term):
    if term["kind"] == "known_unknown":
        return ("left", term["first"], term["second"])
    if term["kind"] == "unknown_known":
        return ("right", term["first"], term["second"])
    if term["kind"] == "unknown":
        return ("identity", term["first"])
    return ("constant", term["first"])


def _canonical_terms(equation):
    result = {}
    for term in equation["terms"]:
        key = _direct_term_key(term)
        result[key] = result.get(key, Fraction(0)) + term["coefficient"]
    return {key: value for key, value in result.items() if value != 0}


def _unknown_name(key):
    if key[0] == "left":
        return key[2]
    if key[0] == "right":
        return key[1]
    if key[0] == "identity":
        return key[1]
    return None


def _scale_matrix(matrix, scalar):
    return [[scalar * value for value in row] for row in matrix]


def _combine_matrices(first, second, first_scale, second_scale):
    return [
        [first[row][col] * first_scale + second[row][col] * second_scale for col in range(len(first[0]))]
        for row in range(len(first))
    ]


def _multiply_matrices(first, second):
    return [
        [
            sum((first[row][k] * second[k][col] for k in range(len(second))), Fraction(0))
            for col in range(len(second[0]))
        ]
        for row in range(len(first))
    ]


def _is_scalar_identity(matrix):
    if not matrix or len(matrix) != len(matrix[0]):
        return False
    scalar = matrix[0][0]
    return scalar != 0 and all(
        matrix[row][col] == (scalar if row == col else 0)
        for row in range(len(matrix))
        for col in range(len(matrix))
    )


def _invert_matrix_silently(matrix):
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        return None
    augmented = [
        matrix[row][:] + [Fraction(1) if row == col else Fraction(0) for col in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column] != 0),
            None,
        )
        if pivot is None:
            return None
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column or augmented[row][column] == 0:
                continue
            multiplier = augmented[row][column]
            augmented[row] = [
                augmented[row][index] - multiplier * augmented[column][index]
                for index in range(2 * size)
            ]
    return [row[size:] for row in augmented]


def _direct_term_latex(key, coefficient):
    if key[0] == "left":
        body = _name_latex(key[1]) + r"\cdot " + _name_latex(key[2])
    elif key[0] == "right":
        body = _name_latex(key[1]) + r"\cdot " + _name_latex(key[2])
    elif key[0] == "identity":
        body = _name_latex(key[1])
    else:
        body = _name_latex(key[1])
    if abs(coefficient) != 1:
        body = latex(abs(coefficient)) + r"\," + body
    return ("-" if coefficient < 0 else "") + body


def _matrix_latex(matrix):
    rows = [r" & ".join(latex(value) for value in row) for row in matrix]
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"


def _combination_latex(first, second, equations):
    first_left = rf"\left({latex(first)}\right)\cdot\left({_equation_lhs_latex(equations[0])}\right)"
    first_right = rf"\left({latex(first)}\right)\cdot {_equation_rhs_latex(equations[0])}"
    second_left = rf"\left({latex(abs(second))}\right)\cdot\left({_equation_lhs_latex(equations[1])}\right)"
    second_right = rf"\left({latex(abs(second))}\right)\cdot {_equation_rhs_latex(equations[1])}"
    if second < 0:
        return first_left + " - " + second_left + " = " + first_right + " - " + second_right
    return first_left + " + " + second_left + " = " + first_right + " + " + second_right


def _direct_term_solution(key, coefficient, rhs, known):
    if coefficient == 0:
        return None
    if key[0] == "identity":
        return {
            "value": _scale_matrix(rhs, Fraction(1, 1) / coefficient),
            "inverse": None,
            "matrix": None,
        }
    known_name = key[1] if key[0] == "left" else key[2]
    coefficient_matrix = _scale_matrix(known[known_name], coefficient)
    if key[0] == "right":
        coefficient_matrix = transpose(coefficient_matrix)
    if not coefficient_matrix or len(coefficient_matrix) != len(coefficient_matrix[0]):
        return None
    if _is_scalar_identity(coefficient_matrix):
        return {
            "value": _scale_matrix(rhs, Fraction(1, 1) / coefficient_matrix[0][0]),
            "inverse": None,
            "matrix": coefficient_matrix,
        }
    inverse = _invert_matrix_silently(coefficient_matrix)
    if inverse is None:
        return None
    if key[0] == "left":
        value = _multiply_matrices(inverse, rhs)
    else:
        value = transpose(_multiply_matrices(inverse, transpose(rhs)))
    return {"value": value, "inverse": inverse, "matrix": coefficient_matrix}


def _direct_addition_plan(prepared, names, known):
    if len(prepared) != 2 or len(names) != 2:
        return None
    equations = [_canonical_terms(equation) for equation in prepared]
    entries = {}
    for name in names:
        rows = []
        for equation in equations:
            matches = [
                (key, coefficient)
                for key, coefficient in equation.items()
                if _unknown_name(key) == name
            ]
            if len(matches) != 1 or matches[0][1] == 0:
                return None
            rows.append(matches[0])
        if rows[0][0] != rows[1][0]:
            return None
        entries[name] = rows
    plan = []
    plan_names = [name for name in ("X", "Y") if name in names]
    plan_names.extend(name for name in names if name not in ("X", "Y"))
    for target in plan_names:
        other = names[1] if target == names[0] else names[0]
        target_key, target_first = entries[target][0]
        target_key_second, target_second = entries[target][1]
        _, other_first = entries[other][0]
        other_key_second, other_second = entries[other][1]
        if target_key_second != target_key or other_key_second != entries[other][0][0]:
            return None
        first_multiplier = other_second
        second_multiplier = -other_first
        if first_multiplier < 0:
            first_multiplier = -first_multiplier
            second_multiplier = -second_multiplier
        combined = {}
        for key, coefficient in equations[0].items():
            combined[key] = combined.get(key, Fraction(0)) + first_multiplier * coefficient
        for key, coefficient in equations[1].items():
            combined[key] = combined.get(key, Fraction(0)) + second_multiplier * coefficient
        combined = {key: coefficient for key, coefficient in combined.items() if coefficient != 0}
        unknown_terms = [
            (key, coefficient)
            for key, coefficient in combined.items()
            if _unknown_name(key) is not None
        ]
        if len(unknown_terms) != 1 or unknown_terms[0][0] != target_key:
            return None
        rhs = _combine_matrices(
            prepared[0]["rhs"],
            prepared[1]["rhs"],
            first_multiplier,
            second_multiplier,
        )
        for key, coefficient in combined.items():
            if key[0] == "constant":
                rhs = _combine_matrices(rhs, known[key[1]], Fraction(1), -coefficient)
        coefficient = unknown_terms[0][1]
        solution = _direct_term_solution(target_key, coefficient, rhs, known)
        if solution is None:
            return None
        plan.append({
            "name": target,
            "other": other,
            "first_multiplier": first_multiplier,
            "second_multiplier": second_multiplier,
            "key": target_key,
            "coefficient": coefficient,
            "rhs": rhs,
            "solution": solution,
        })
    return plan


def _write_direct_addition(plan, equations, w):
    w.h("Метод алгебраического сложения", 2)
    for item in plan:
        target = item["name"]
        other = item["other"]
        first = item["first_multiplier"]
        second = item["second_multiplier"]
        w.h(f"Сокращаем {other} и находим {target}", 2)
        w.p(f"Умножаем первое уравнение на {latex(first)}.")
        w.l(_scaled_equation_latex(equations[0], first))
        w.p(f"Умножаем второе уравнение на {latex(second)}.")
        w.l(_scaled_equation_latex(equations[1], second))
        w.p(f"Складываем полученные уравнения и сокращаем {other}.")
        w.l(_combination_latex(first, second, equations))
        w.l(
            f"{_direct_term_latex(item['key'], item['coefficient'])} = "
            f"{_matrix_latex(item['rhs'])}"
        )
        solution = item["solution"]
        if solution["inverse"] is None:
            w.p(f"Сокращённое уравнение сразу даёт {target} = найденную матрицу.")
        else:
            known_name = item["key"][1] if item["key"][0] == "left" else item["key"][2]
            coefficient_title = _name_latex(known_name) + (r"^{T}" if item["key"][0] == "right" else "")
            w.m(
                "Коэффициент",
                solution["matrix"],
                caption=f"Матрица коэффициента для {coefficient_title}",
            )
            inverse_title = coefficient_title + r"^{-1}"
            w.m(inverse_title, solution["inverse"], caption="Обратная матрица коэффициента")
            w.p(f"Умножаем на {inverse_title} и получаем {target} = найденную матрицу.")
        w.l(
            rf"\boxed{{{target} = {_matrix_latex(solution['value'])}}}",
            ans=True,
        )
        w.m(target, solution["value"], caption=f"Матрица {target}", ans=True)


def _equation_lhs_latex(equation, scalar=None):
    parts = []
    for term in equation["terms"]:
        coefficient = term["coefficient"] if scalar is None else term["coefficient"] * scalar
        if term["kind"] == "known_unknown":
            body = _name_latex(term["first"]) + r"\cdot " + _name_latex(term["second"])
        elif term["kind"] == "unknown_known":
            body = _name_latex(term["first"]) + r"\cdot " + _name_latex(term["second"])
        else:
            body = _name_latex(term["first"])
        if abs(coefficient) != 1:
            body = latex(abs(coefficient)) + r"\," + body
        if not parts:
            parts.append(("-" if coefficient < 0 else "") + body)
        else:
            parts.append((" + " if coefficient > 0 else " - ") + body)
    return "".join(parts) if parts else "0"


def _equation_rhs_latex(equation, scalar=None):
    rhs = _name_latex(equation["rhs_name"])
    if scalar is not None and abs(scalar) != 1:
        rhs = latex(abs(scalar)) + r"\," + rhs
        if scalar < 0:
            rhs = "-" + rhs
    return rhs


def _equation_latex(equation):
    return _equation_lhs_latex(equation) + " = " + _equation_rhs_latex(equation)


def _scaled_equation_latex(equation, scalar):
    return (
        _equation_lhs_latex(equation, scalar)
        + " = "
        + _equation_rhs_latex(equation, scalar)
    )


def _system_latex(equations):
    lines = [_equation_latex(equation) for equation in equations]
    return r"\left\{\begin{aligned}" + r" \\ ".join(lines) + r"\end{aligned}\right."



def _parameter(index, group=None):
    if group is None:
        return "t_{" + str(index + 1) + "}"
    return "t_{" + str(group + 1) + "," + str(index + 1) + "}"


def _expression(constant, coefficients, parameters):
    parts = []
    if constant != 0:
        parts.append(latex(constant))
    for coefficient, parameter in zip(coefficients, parameters):
        if coefficient == 0:
            continue
        absolute = abs(coefficient)
        if absolute == 1:
            term = parameter
        else:
            term = latex(absolute) + r"\," + parameter
        if not parts:
            parts.append(("-" if coefficient < 0 else "") + term)
        else:
            parts.append((" + " if coefficient > 0 else " - ") + term)
    return "".join(parts) if parts else "0"


def _expressions(matrix, pivot_rows, free, n, q, layout, grouped):
    result = [["" for _ in range(0)] for _ in range(0)]
    for group in range(q):
        values = []
        for index in range(n):
            if index in pivot_rows:
                row = pivot_rows[index]
                constant = matrix[row][n + group]
                coefficients = [matrix[row][free_col] for free_col in free]
                parameters = [_parameter(free_col, group if grouped else None) for free_col in free]
            else:
                constant = Fraction(0)
                coefficients = [Fraction(1) if free_col == index else Fraction(0) for free_col in free]
                parameters = [_parameter(free_col, group if grouped else None) for free_col in free]
            values.append(_expression(constant, coefficients, parameters))
        result.append(values)
    return result


def _unique_values(matrix, pivot_rows, n, q):
    values = []
    for group in range(q):
        column = []
        for index in range(n):
            if index not in pivot_rows:
                column.append(None)
            else:
                column.append(matrix[pivot_rows[index]][n + group])
        values.append(column)
    return values


def _layout_matrix(rows, cols):
    return [(row, col) for row in range(rows) for col in range(cols)]


def _layout_named(shapes, names):
    layout = []
    for name in names:
        rows, cols = shapes[name]
        for row in range(rows):
            for col in range(cols):
                layout.append((name, row, col))
    return layout


def _matrix_from_values(values, rows, cols, offset=0):
    return [
        [values[offset + row * cols + col] for col in range(cols)]
        for row in range(rows)
    ]


def _raw_matrix_from_values(values, rows, cols, offset=0):
    return [
        [values[offset + row * cols + col] for col in range(cols)]
        for row in range(rows)
    ]


def _matrix_from_columns(columns, rows, cols):
    return [[columns[col][row] for col in range(cols)] for row in range(rows)]


def solve_axb(A, B, w, unknown_name="X"):
    if not is_identifier(unknown_name):
        raise ValueError("Имя неизвестной матрицы должно начинаться с буквы.")
    if not A or not B:
        raise ValueError("Для матричного уравнения нужны матрицы A и B.")
    a_rows, a_cols = _shape(A)
    b_rows, b_cols = _shape(B)
    if a_cols < 1 or b_cols < 1:
        raise ValueError("Матрицы A и B должны иметь хотя бы один столбец.")
    if a_rows != b_rows:
        raise ValueError(
            f"Число строк A ({a_rows}) не равно числу строк B ({b_rows})."
        )
    w.h(f"Решение матричного уравнения A·{unknown_name} = B", 1)
    w.m("A", A, caption=f"Коэффициентная матрица A ({a_rows}×{a_cols})")
    w.m("B", B, caption=f"Правая матрица B ({a_rows}×{b_cols})")
    w.p(
        f"Ищем матрицу {unknown_name} размера {a_cols}×{b_cols} из уравнения A·{unknown_name} = B. "
        "Приведём расширенную матрицу к ступенчатому виду."
    )
    C = [row[:] for row in A]
    D = [row[:] for row in B]
    matrix, pivots, pivot_rows, free, consistent = _rref_with_steps(C, D, w, "[A | B]")
    rank_augmented = sum(any(value != 0 for value in row) for row in matrix)
    w.p(f"Ранг A = {len(pivots)}, ранг расширенной матрицы = {rank_augmented}.")
    if not consistent:
        w.p("Система несовместна — матричное уравнение не имеет решений.", ans=True)
        w.l(r"\boxed{\text{Матричное уравнение не имеет решений}}", ans=True)
        return None
    layout = _layout_matrix(a_cols, b_cols)
    if not free:
        values = _unique_values(matrix, pivot_rows, a_cols, b_cols)
        X = _matrix_from_columns(values, a_cols, b_cols)
        w.h("Решение", 1)
        w.m(unknown_name, X, caption=f"Матрица {unknown_name} ({a_cols}×{b_cols})", ans=True)
        return X
    expressions = _expressions(matrix, pivot_rows, free, a_cols, b_cols, layout, True)
    X = _matrix_from_columns(expressions, a_cols, b_cols)
    w.p(
        "Есть свободные элементы X — матричное уравнение имеет бесконечно много решений. "
        "Параметры t независимы для разных правых столбцов."
    )
    w.h("Решение (в параметрическом виде)", 1)
    w.mraw(unknown_name, X, caption=f"Матрица {unknown_name} (t — произвольные параметры)", ans=True)
    return X


def solve_system(req, w):
    known = _parse_known(req.get("known", req.get("matrices")))
    shapes = _parse_unknowns(req.get("unknowns"))
    equations = parse_equations(req.get("equations"))
    prepared, names = _prepare_equations(equations, known, shapes)
    direct_plan = _direct_addition_plan(prepared, names, known)
    w.h("Решение системы матричных уравнений", 1)
    if direct_plan is None:
        w.p(
            "Каждое уравнение разворачиваем по элементам: неизвестные матрицы "
            "становятся скалярными переменными, а коэффициенты образуют общую систему."
        )
    else:
        w.p(
            "Сначала складываем и вычитаем исходные матричные уравнения: "
            "сокращаем Y, находим X, затем сокращаем X и находим Y."
        )
    for name, matrix in known.items():
        w.m(name, matrix, caption=f"Известная матрица {name} ({len(matrix)}×{len(matrix[0])})")
    for name in names:
        rows, cols = shapes[name]
        w.p(f"Неизвестная матрица {name}: {rows}×{cols}.")
    w.l(_system_latex(prepared))
    if direct_plan is not None:
        _write_direct_addition(direct_plan, prepared, w)
        return {
            item["name"]: item["solution"]["value"]
            for item in direct_plan
        }
    coefficient_rows, values = _build_scalar_system(prepared, shapes, names, known)
    C = coefficient_rows
    D = [[value] for value in values]
    w.h("Сведение к скалярной системе", 2)
    w.m("[C | d]", _augmented_rows(C, D, True), caption="Каждое уравнение соответствует одному элементу")
    w.h("Метод алгебраического сложения", 2)
    w.p(
        "Подбираем линейные комбинации уравнений, чтобы сократить очередной неизвестный, "
        "а затем прибавляем или вычитаем полученные уравнения."
    )
    matrix, pivots, pivot_rows, free, consistent = _rref_by_addition(C, D, w, "[C | d]")
    rank_augmented = sum(any(value != 0 for value in row) for row in matrix)
    w.p(f"Ранг C = {len(pivots)}, ранг расширенной матрицы = {rank_augmented}.")
    if not consistent:
        w.p("Система несовместна — решений нет.", ans=True)
        w.l(r"\boxed{\text{Система несовместна — решений нет}}", ans=True)
        return None
    layout = _layout_named(shapes, names)
    if not free:
        values_flat = _unique_values(matrix, pivot_rows, len(layout), 1)[0]
        offset = 0
        w.h("Решение", 1)
        for name in names:
            rows, cols = shapes[name]
            current = _matrix_from_values(values_flat, rows, cols, offset)
            w.m(name, current, caption=f"Матрица {name} ({rows}×{cols})", ans=True)
            offset += rows * cols
        return {name: _matrix_from_values(values_flat, shapes[name][0], shapes[name][1], sum(shapes[previous][0] * shapes[previous][1] for previous in names[:names.index(name)])) for name in names}
    expressions_flat = _expressions(matrix, pivot_rows, free, len(layout), 1, layout, False)[0]
    w.p(
        "Есть свободные переменные — система имеет бесконечно много решений. "
        "Параметры t независимы для разных элементов неизвестных матриц."
    )
    w.h("Решение (в параметрическом виде)", 1)
    result = {}
    offset = 0
    for name in names:
        rows, cols = shapes[name]
        current = _raw_matrix_from_values(expressions_flat, rows, cols, offset)
        w.mraw(name, current, caption=f"Матрица {name} (t — произвольные параметры)", ans=True)
        result[name] = current
        offset += rows * cols
    return result
