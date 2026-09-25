import re


def _build_matrix(data):
    sepAt = -1
    for row in data:
        cnt = 0
        for c in row:
            if isinstance(c, dict) and c.get("sep"):
                sepAt = cnt
                break
            cnt += 1
        if sepAt >= 0:
            break

    rows = []
    for row in data:
        cells = [c for c in row if not (isinstance(c, dict) and c.get("sep"))]
        rows.append(" & ".join(str(c) for c in cells))

    if sepAt >= 0:
        total = len([c for c in data[0] if not (isinstance(c, dict) and c.get("sep"))])
        spec = "c" * sepAt + "|" + "c" * max(total - sepAt, 0)
        return (
            "\\left[\\begin{array}{" + spec + "}"
            + " \\\\ ".join(rows)
            + "\\end{array}\\right]"
        )
    return "\\begin{pmatrix}" + " \\\\ ".join(rows) + "\\end{pmatrix}"


def _is_math(text):
    return any(ch in text for ch in ("\\", "_", "{", "}", "^"))


def _strip_text(line):
    return re.sub(r"\\text\{[^{}]*\}", "", line).strip()


def build_latex(steps):
    doc = []
    doc.append("\\documentclass[12pt]{article}")
    doc.append("\\usepackage[utf8]{inputenc}")
    doc.append("\\usepackage[T2A]{fontenc}")
    doc.append("\\usepackage[russian]{babel}")
    doc.append("\\usepackage[a4paper,margin=2cm]{geometry}")
    doc.append("\\usepackage{amsmath,amssymb}")
    doc.append("\\pagestyle{plain}")
    doc.append("")
    doc.append("\\begin{document}")

    for s in steps:
        t = s.get("t")
        if t in ("h", "p"):
            continue
        if t == "l":
            doc.append("\\[" + _strip_text(s.get("s", "")) + "\\]")
        elif t == "m":
            data = s.get("data") or []
            title = s.get("title")
            mx = _build_matrix(data)
            if title and _is_math(str(title)):
                doc.append(f"\\[{title} = {mx}\\]")
            else:
                doc.append("\\[" + mx + "\\]")

    doc.append("\\end{document}")
    return "\n".join(doc)