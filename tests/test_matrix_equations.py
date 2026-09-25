import re
import unittest
from fractions import Fraction
from pathlib import Path

from solver.api import calc
from solver.equation_parser import parse_equations, parse_matrix_equation

ROOT = Path(__file__).resolve().parent.parent


class EquationParserTests(unittest.TestCase):
    def test_parses_linear_matrix_expression(self):
        parsed = parse_equations("2*A*X - 1/2*B·Y = C")
        self.assertEqual(parsed[0]["rhs"], "C")
        self.assertEqual(parsed[0]["terms"][0]["coefficient"], Fraction(2))
        self.assertEqual(parsed[0]["terms"][0]["first"], "A")
        self.assertEqual(parsed[0]["terms"][0]["second"], "X")
        self.assertEqual(parsed[0]["terms"][1]["coefficient"], Fraction(-1, 2))

    def test_accepts_several_lines_and_structured_terms(self):
        parsed = parse_equations([
            "A*X = B\nA*Y = C",
            {"terms": [{"known": "A", "unknown": "X"}], "rhs": "B"},
        ])
        self.assertEqual(len(parsed), 3)
        self.assertEqual(parsed[2]["terms"][0]["second"], "X")

    def test_parses_inline_matrix_equation(self):
        parsed = parse_matrix_equation("[[2,0],[0,3]]*Y=[[4,2],[0,9]]")
        self.assertEqual(parsed["unknown"], "Y")
        self.assertEqual(parsed["A"], [[Fraction(2), Fraction(0)], [Fraction(0), Fraction(3)]])
        self.assertEqual(parsed["B"], [[Fraction(4), Fraction(2)], [Fraction(0), Fraction(9)]])

    def test_accepts_every_transpose_spelling(self):
        for text in ("AT*X = B", "A^T*X = B", "A'*X = B", "A\u1d40*X = B", "A\u1d57*X = B"):
            parsed = parse_equations(text, {"A", "B"})[0]
            self.assertIn(r"A^{T}\cdot X = B", parsed["latex"], text)
            self.assertEqual(parsed["terms"][0]["first"], "A", text)
            self.assertTrue(parsed["terms"][0]["first_transposed"], text)
            self.assertEqual(parsed["terms"][0]["second"], "X", text)

    def test_transposes_parenthesised_expression(self):
        parsed = parse_equations("(2A - 3C)T*X = B", {"A", "B", "C"})[0]
        self.assertIn(r"\left(2\,A - 3\,C\right)^{T}\cdot X = B", parsed["latex"])
        factor = parsed["terms"][0]["factors"][0]
        self.assertEqual(factor["kind"], "expr")
        self.assertTrue(factor["transposed"])

    def test_transposes_right_hand_side(self):
        parsed = parse_equations("A*X = (B + C)^T", {"A", "B", "C"})[0]
        self.assertIn(r"X = \left(B + C\right)^{T}", parsed["latex"])
        self.assertTrue(parsed["rhs_terms"][0]["factors"][0]["transposed"])
        parsed = parse_equations("A*X = C^T", {"A", "C"})[0]
        self.assertEqual(parsed["rhs"], "C")
        self.assertTrue(parsed["rhs_transposed"])

    def test_flips_order_when_transposing_a_product(self):
        parsed = parse_equations("(A*B)^T*X = C", {"A", "B", "C"})[0]
        self.assertIn(r"B^{T}\cdot A^{T}\cdot X = C", parsed["latex"])

    def test_keeps_declared_matrix_name_ending_with_t(self):
        parsed = parse_equations("CT*X = B", {"CT", "B"})[0]
        self.assertEqual(parsed["terms"][0]["first"], "CT")
        self.assertFalse(parsed["terms"][0]["first_transposed"])
        self.assertIn(r"CT\cdot X = B", parsed["latex"])

    def test_rejects_transpose_of_unknown(self):
        result = calc({
            "op": "matrixsys",
            "known": {"A": [[1, 0], [0, 1]], "B": [[1, 0], [0, 1]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "A*X^T = B",
        })
        self.assertIn("error", result)
        self.assertIn("не поддерживается", result["error"])

    def test_rejects_unknown_inside_parentheses(self):
        result = calc({
            "op": "matrixsys",
            "known": {"A": [[1, 0], [0, 1]], "B": [[1, 0], [0, 1]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "(A*X + B) = B",
        })
        self.assertIn("error", result)
        self.assertIn("только известные", result["error"])

    def test_rejects_missing_transpose_marker(self):
        with self.assertRaises(ValueError):
            parse_equations("A^X = B", {"A", "B"})


class MatrixEquationTests(unittest.TestCase):
    @staticmethod
    def last_matrix(result, title):
        for step in reversed(result["steps"]):
            if step.get("t") == "m" and step.get("title") == title and step.get("ans"):
                return step["data"]
        raise AssertionError(f"Не найден ответ {title}: {result}")

    def test_solves_axb_with_two_right_hand_sides(self):
        result = calc({
            "op": "matrixeq",
            "A": [[2, 0], [0, 3]],
            "B": [[2, 4], [3, 6]],
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "2"], ["1", "2"]])

    def test_solves_inline_matrix_equation(self):
        result = calc({
            "op": "matrixeq",
            "equation": "[[2,0],[0,3]]*Y=[[4,2],[0,9]]",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "Y"), [["2", "1"], ["0", "3"]])

    def test_solves_rectangular_axb(self):
        result = calc({
            "op": "matrixeq",
            "A": [[1, 0], [0, 1], [1, 1]],
            "B": [[1, 2], [3, 4], [4, 6]],
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "2"], ["3", "4"]])

    def test_reports_inconsistent_axb(self):
        result = calc({
            "op": "matrixeq",
            "A": [[1, 0], [0, 0]],
            "B": [[1], [2]],
        })
        self.assertNotIn("error", result)
        self.assertTrue(any("несовместна" in step.get("s", "") for step in result["steps"]))

    def test_returns_parametric_axb(self):
        result = calc({
            "op": "matrixeq",
            "A": [[1, 1, 0]],
            "B": [[2, 3]],
        })
        self.assertNotIn("error", result)
        data = self.last_matrix(result, "X")
        self.assertIn("t_{1,2}", data[0][0])
        self.assertEqual(data[1][0], "t_{1,2}")

    def test_latex_mode(self):
        result = calc({
            "op": "matrixeq",
            "A": [[1]],
            "B": [[2]],
            "mode": "latex",
        })
        self.assertIn("latex", result)
        self.assertIn("X", result["latex"])

    def test_latex_preserves_inconsistent_conclusion(self):
        result = calc({
            "op": "matrixeq",
            "A": [[0]],
            "B": [[1]],
            "mode": "latex",
        })
        self.assertIn(r"\text{Матричное уравнение не имеет решений}", result["latex"])


class MatrixSystemTests(unittest.TestCase):
    @staticmethod
    def last_matrix(result, title):
        for step in reversed(result["steps"]):
            if step.get("t") == "m" and step.get("title") == title and step.get("ans"):
                return step["data"]
        raise AssertionError(f"Не найден ответ {title}: {result}")

    def test_solves_two_matrix_unknowns(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 0], [0, 1]],
                "B": [[1, 0], [0, 1]],
                "C": [[3, 1], [2, 4]],
                "D": [[1, 0], [0, 1]],
            },
            "unknowns": [
                {"name": "X", "rows": 2, "cols": 2},
                {"name": "Y", "rows": 2, "cols": 2},
            ],
            "equations": "A*X + B*Y = C\nA*X - B*Y = D",
        })
        self.assertNotIn("error", result)
        self.assertTrue(any(r"\left\{\begin{aligned}" in step.get("s", "") for step in result["steps"]))
        self.assertTrue(any(step.get("s") == "Метод алгебраического сложения" for step in result["steps"]))
        self.assertTrue(
            any(step.get("s") == "Сокращаем Y и находим X" for step in result["steps"])
        )
        self.assertTrue(
            any(step.get("s") == "Сокращаем X и находим Y" for step in result["steps"])
        )
        self.assertTrue(
            any(
                r"\left(1\right)\cdot\left(A\cdot X + B\cdot Y\right)" in step.get("s", "")
                for step in result["steps"]
            )
        )
        self.assertTrue(
            any(
                r"\left(1\right)\cdot\left(A\cdot X - B\cdot Y\right)" in step.get("s", "")
                for step in result["steps"]
            )
        )
        self.assertTrue(any(r"2\,A\cdot X =" in step.get("s", "") for step in result["steps"]))
        self.assertTrue(any(r"\boxed{X =" in step.get("s", "") for step in result["steps"]))
        self.assertTrue(any(r"\boxed{Y =" in step.get("s", "") for step in result["steps"]))
        self.assertFalse(any(step.get("s") == "Сведение к скалярной системе" for step in result["steps"]))
        self.assertFalse(any("E_" in step.get("s", "") for step in result["steps"]))
        self.assertEqual(self.last_matrix(result, "X"), [["2", r"\frac{1}{2}"], ["1", r"\frac{5}{2}"]])
        self.assertEqual(self.last_matrix(result, "Y"), [["1", r"\frac{1}{2}"], ["1", r"\frac{3}{2}"]])

    def test_matches_coefficients_before_addition(self):
        result = calc({
            "op": "matrixsys",
            "known": {"C": [[1]], "D": [[3]]},
            "unknowns": [
                {"name": "X", "rows": 1, "cols": 1},
                {"name": "Y", "rows": 1, "cols": 1},
            ],
            "equations": "2*X + 3*Y = C\n4*X + 5*Y = D",
        })
        self.assertNotIn("error", result)
        self.assertTrue(any(step.get("s") == "Сокращаем Y и находим X" for step in result["steps"]))
        self.assertTrue(
            any(
                r"\left(5\right)\cdot\left(2\,X + 3\,Y\right) - \left(3\right)\cdot\left(4\,X + 5\,Y\right)"
                in step.get("s", "")
                for step in result["steps"]
            )
        )
        self.assertTrue(any(r"10\,X + 15\,Y = 5\,C" in step.get("s", "") for step in result["steps"]))
        self.assertFalse(any(step.get("s") == "Сведение к скалярной системе" for step in result["steps"]))
        self.assertFalse(any("E_" in step.get("s", "") for step in result["steps"]))
        self.assertEqual(self.last_matrix(result, "X"), [["2"]])
        self.assertEqual(self.last_matrix(result, "Y"), [["-1"]])

    def test_solves_direct_addition_with_matrix_coefficients(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 1], [0, 1]],
                "B": [[1, 0], [1, 1]],
                "C": [[3, 6], [5, 5]],
                "D": [[1, 2], [-5, 1]],
            },
            "unknowns": [
                {"name": "X", "rows": 2, "cols": 2},
                {"name": "Y", "rows": 2, "cols": 2},
            ],
            "equations": "A*X + B*Y = C\nA*X - B*Y = D",
        })
        self.assertNotIn("error", result)
        self.assertFalse(any(step.get("s") == "Сведение к скалярной системе" for step in result["steps"]))
        self.assertFalse(any("E_" in step.get("s", "") for step in result["steps"]))
        self.assertTrue(any("A^{-1}" in step.get("s", "") for step in result["steps"]))
        self.assertEqual(self.last_matrix(result, "X"), [["2", "1"], ["0", "3"]])
        self.assertEqual(self.last_matrix(result, "Y"), [["1", "2"], ["4", "0"]])

    def test_supports_unknown_on_left(self):
        result = calc({
            "op": "matrixsys",
            "known": {"A": [[1, 0], [0, 1]], "B": [[1, 2], [3, 4]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "X*A = B",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "2"], ["3", "4"]])

    def test_supports_rectangular_unknown_on_left(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 0], [0, 1], [1, 1]],
                "E": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "B": [[4, 5], [10, 11]],
                "F": [[1, 2, 3], [4, 5, 6]],
            },
            "unknowns": [{"name": "X"}],
            "equations": "X*A = B\nX*E = F",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "2", "3"], ["4", "5", "6"]])

    def test_accepts_known_matrix_list_from_ui(self):
        result = calc({
            "op": "matrixsys",
            "known": [{"name": "A", "data": [[2]]}],
            "unknowns": [{"name": "X", "rows": 1, "cols": 1}],
            "equations": "A*X = A",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1"]])

    def test_solves_equation_with_transposed_matrices(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 1], [0, 1]],
                "C": [[1, 0], [0, 2]],
            },
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "5X + 3A - 2C^T = (2A - 3C)T",
        })
        self.assertNotIn("error", result)
        self.assertTrue(
            any(r"5\,X + 3\,A - 2\,C^{T} = \left(2\,A - 3\,C\right)^{T}" in step.get("s", "")
                for step in result["steps"])
        )
        self.assertFalse(any(r"\mathrm{" in step.get("s", "") for step in result["steps"]))
        self.assertEqual(
            self.last_matrix(result, "X"),
            [[r"\frac{-2}{5}", r"\frac{-3}{5}"], [r"\frac{2}{5}", r"\frac{-3}{5}"]],
        )

    def test_solves_system_with_transposed_coefficients(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 1], [0, 1]],
                "B": [[1, 0], [1, 1]],
                "C": [[3, 6], [5, 5]],
                "D": [[1, 2], [-5, 1]],
            },
            "unknowns": [
                {"name": "X", "rows": 2, "cols": 2},
                {"name": "Y", "rows": 2, "cols": 2},
            ],
            "equations": "A^T*X + B^T*Y = C\nA^T*X - B^T*Y = D",
        })
        self.assertNotIn("error", result)
        self.assertFalse(any(step.get("s") == "Сведение к скалярной системе" for step in result["steps"]))
        self.assertEqual(self.last_matrix(result, "X"), [["2", "4"], ["-2", "-1"]])
        self.assertEqual(self.last_matrix(result, "Y"), [["-4", "0"], ["5", "2"]])

    def test_solves_unknown_between_two_transposed_matrices(self):
        result = calc({
            "op": "matrixsys",
            "known": {
                "A": [[1, 2], [3, 4]],
                "B": [[2, 0], [0, 1]],
                "C": [[2, 3], [4, 4]],
            },
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "A^T*X*B^T = C",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "0"], ["0", "1"]])

    def test_matrixeq_op_accepts_known_and_equations(self):
        result = calc({
            "op": "matrixeq",
            "known": {"A": [[2, 0], [0, 3]], "B": [[2, 4], [3, 6]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "A*X = B",
        })
        self.assertNotIn("error", result)
        self.assertEqual(self.last_matrix(result, "X"), [["1", "2"], ["1", "2"]])

    def test_matrixeq_op_reports_empty_payload(self):
        result = calc({"op": "matrixeq"})
        self.assertIn("error", result)
        self.assertIn("известные матрицы", result["error"])

    def test_rejects_name_collision(self):
        result = calc({
            "op": "matrixsys",
            "known": {"X": [[1]], "B": [[1]]},
            "unknowns": [{"name": "X", "rows": 1, "cols": 1}],
            "equations": "X = B",
        })
        self.assertIn("error", result)
        self.assertIn("известную", result["error"])

    def test_reports_inconsistent_system(self):
        result = calc({
            "op": "matrixsys",
            "known": {"A": [[0, 0], [0, 0]], "B": [[1, 0], [0, 1]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "A*X = B",
        })
        self.assertNotIn("error", result)
        self.assertTrue(any("несовместна" in step.get("s", "") for step in result["steps"]))

    def test_validates_matrix_dimensions(self):
        result = calc({
            "op": "matrixsys",
            "known": {"A": [[1, 0], [0, 1]], "B": [[1, 0], [0, 1]]},
            "unknowns": [{"name": "X", "rows": 2, "cols": 2}],
            "equations": "A*X = B",
        })
        self.assertNotIn("error", result)


class FrontendWiringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.html = (ROOT / "static" / "index.html").read_text(encoding="utf-8")
        cls.js = (ROOT / "static" / "app.js").read_text(encoding="utf-8")

    def test_every_element_id_used_by_script_exists_in_page(self):
        ids = set(re.findall(r'\$\("([A-Za-z][\w-]*)"\)', self.js))
        missing = sorted(name for name in ids if f'id="{name}"' not in self.html)
        self.assertEqual(missing, [])

    def test_matrix_equation_and_system_share_one_editor(self):
        self.assertEqual(self.js.count("matrixeq: { system: true, single: true }"), 1)
        self.assertNotIn("matrixEquation", self.js)
        self.assertNotIn("matrixEquation", self.html)
        self.assertIn('id="sec-matrixsys"', self.html)
        self.assertIn("meta.single", self.js)

    def test_example_uses_transposed_equation(self):
        self.assertIn("5X + 3A - 2C^T = (2A - 3C)T", self.js)
        self.assertIn("(2A - 3C)T", self.html)

    def test_no_operation_needs_the_inline_equation_field(self):
        for op in re.findall(r"^\s*(\w+):\s*\{([^}]*)\},", self.js, re.M):
            self.assertNotIn("equation: true", op[1], op[0])


if __name__ == "__main__":
    unittest.main()
