from . import writer
from .frac import parse_num
from . import ops, det, minors, inverse, slae

MAX_SIZE = 6

SQUARE_OPS = {"det", "detg", "minors", "invadj", "invgauss", "cramer", "gauss"}
NEED_B = {"cramer", "gauss"}


def parse_mat(data):
    if not isinstance(data, list) or len(data) == 0:
        raise ValueError("Матрица не задана.")
    M = []
    for row in data:
        if not isinstance(row, list):
            row = [row]
        M.append([parse_num(v) for v in row])
    if len(M) > MAX_SIZE or any(len(row) > MAX_SIZE for row in M):
        raise ValueError(f"Размер матриц ограничен {MAX_SIZE}×{MAX_SIZE}.")
    if len({len(row) for row in M}) > 1:
        raise ValueError("Все строки матрицы должны быть одинаковой длины.")
    return M


def normalize_b(bs):
    if not isinstance(bs, list):
        raise ValueError("Столбец свободных членов b не задан.")
    out = []
    for x in bs:
        if isinstance(x, list):
            if len(x) != 1:
                raise ValueError("Каждая строка столбца b должна содержать одно значение.")
            out.append(x[0])
        else:
            out.append(x)
    return [parse_num(v) for v in out]


def _check_square(op, A):
    if op in SQUARE_OPS and len(A) != len(A[0]):
        raise ValueError("Для этой операции матрица должна быть квадратной (n×n).")


def calc(req):
    op = req.get("op")
    w = writer.W()
    try:
        A = parse_mat(req["A"]) if req.get("A") is not None else None
        B = parse_mat(req["B"]) if req.get("B") is not None else None

        if op == "add":
            if A is None or B is None:
                raise ValueError("Нужны обе матрицы A и B.")
            ops.add_sub(A, B, w, "+")
        elif op == "sub":
            if A is None or B is None:
                raise ValueError("Нужны обе матрицы A и B.")
            ops.add_sub(A, B, w, "-")
        elif op == "scalar":
            if A is None:
                raise ValueError("Нужна матрица A.")
            k = parse_num(req.get("k", "1"))
            ops.scalar_mul(k, A, w)
        elif op == "mul":
            if A is None or B is None:
                raise ValueError("Нужны обе матрицы A и B.")
            if not A or not B or len(A[0]) != len(B):
                raise ValueError(
                    f"Матрицы несогласованы: число столбцов A ({len(A[0])}) "
                    f"не равно числу строк B ({len(B)}).")
            ops.mat_mul(A, B, w)
        elif op == "transpose":
            if A is None:
                raise ValueError("Нужна матрица A.")
            ops.transpose_op(A, w)
        elif op == "rank":
            if A is None:
                raise ValueError("Нужна матрица A.")
            ops.rank(A, w)
        elif op in ("det", "detg", "minors", "invadj", "invgauss"):
            if A is None:
                raise ValueError("Нужна матрица A.")
            _check_square(op, A)
            if op == "det":
                exp = req.get("expand") or {"kind": "row", "idx": 0}
                kind = exp.get("kind", "row")
                idx = int(exp.get("idx", 0))
                if kind == "row":
                    if not (0 <= idx < len(A)):
                        raise ValueError("Выбранной строки не существует.")
                    det.laplace(A, w, expand=("row", idx))
                else:
                    if not (0 <= idx < len(A[0])):
                        raise ValueError("Выбранного столбца не существует.")
                    det.laplace(A, w, expand=("col", idx))
            elif op == "detg":
                det.det_gauss(A, w)
            elif op == "minors":
                minors.minors_cofactors(A, w)
            elif op == "invadj":
                inverse.inverse_adjugate(A, w)
            elif op == "invgauss":
                inverse.inverse_gauss(A, w)
        elif op in NEED_B:
            if A is None:
                raise ValueError("Нужна матрица A.")
            b = normalize_b(req.get("b"))
            if len(b) != len(A):
                raise ValueError(f"Число свободных членов ({len(b)}) не равно числу уравнений ({len(A)}).")
            if op == "cramer":
                _check_square(op, A)
                slae.cramer(A, b, w)
            else:
                slae.gauss(A, b, w)
        else:
            raise ValueError("Неизвестная операция: " + str(op))
    except Exception as e:
        return {"error": str(e)}

    steps = w.steps
    if not steps:
        steps = [{"t": "p", "s": "Ничего не вычислено."}]
    return {"steps": steps}