import unittest
from fractions import Fraction

from solver.api import calc
from solver.equation_parser import parse_equations, parse_matrix_equation


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


if __name__ == "__main__":
    unittest.main()
