def _cell(v):
    from .frac import latex

    return latex(v)


def _rows(data):
    out = []
    for row in data:
        r = []
        for v in row:
            if isinstance(v, dict) and v.get("sep"):
                r.append({"sep": True})
            else:
                r.append(_cell(v))
        out.append(r)
    return out


class W:
    def __init__(self):
        self.steps = []

    def h(self, s, lvl=1):
        self.steps.append({"t": "h", "s": s, "l": lvl})

    def p(self, s, ans=False):
        self.steps.append({"t": "p", "s": s, "ans": ans})

    def l(self, s, ans=False):
        self.steps.append({"t": "l", "s": s, "ans": ans})

    def m(self, title, data, caption=None, ans=False):
        self.steps.append({
            "t": "m",
            "title": title,
            "data": _rows(data),
            "caption": caption,
            "ans": ans,
        })

    def mraw(self, title, rows, caption=None, ans=False):
        self.steps.append({
            "t": "m",
            "title": title,
            "data": rows,
            "caption": caption,
            "ans": ans,
        })