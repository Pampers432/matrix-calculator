from .frac import cell


def rows_data(data):
    out = []
    for row in data:
        r = []
        for v in row:
            if isinstance(v, dict) and v.get("sep"):
                r.append({"sep": True})
            else:
                r.append(cell(v))
        out.append(r)
    return out


class W:
    def __init__(self):
        self.steps = []

    def h(self, s, lvl=1):
        self.steps.append({"t": "h", "s": s, "l": lvl})

    def p(self, s, ans=False):
        self.steps.append({"t": "p", "s": s, "ans": ans})

    def m(self, title, data, caption=None, ans=False):
        self.steps.append({
            "t": "m",
            "title": title,
            "data": rows_data(data),
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