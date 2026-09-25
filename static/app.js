"use strict";

const $ = (id) => document.getElementById(id);

const OPS = {
  add:      { b: true },
  sub:      { b: true },
  scalar:   { k: true },
  mul:      { b: true },
  transpose: {},
  det:      { expand: true },
  detg:     {},
  minors:   {},
  rank:     {},
  invadj:   {},
  invgauss: {},
  cramer:   { rhs: true },
  gauss:    { rhs: true },
};

const opSel = $("op");
const rowsA = $("rowsA"), colsA = $("colsA");
const rowsB = $("rowsB"), colsB = $("colsB");
const gridA = $("gridA"), gridB = $("gridB"), gridb = $("gridb");

const EXAMPLE_A = [[1, 2, 3], [2, 1, 3], [3, 2, 1]];
const EXAMPLE_B = [[1, 0, 2], [0, 1, 3], [2, 3, 1]];
const EXAMPLE_b = [[6], [6], [6]];

function buildGrid(container, rows, cols, data) {
  container.innerHTML = "";
  const table = document.createElement("table");
  const cls = rows * cols > 16 ? "small" : (rows * cols > 9 ? "mid" : "");
  for (let r = 0; r < rows; r++) {
    const tr = document.createElement("tr");
    for (let c = 0; c < cols; c++) {
      const td = document.createElement("td");
      const inp = document.createElement("input");
      inp.type = "text";
      inp.className = cls;
      inp.autocomplete = "off";
      inp.spellcheck = false;
      const v = data && data[r] && data[r][c] !== undefined ? data[r][c] : "0";
      inp.value = String(v);
      td.appendChild(inp);
      tr.appendChild(td);
    }
    table.appendChild(tr);
  }
  container.appendChild(table);
}

function readGrid(container) {
  const out = [];
  container.querySelectorAll("tr").forEach((tr) => {
    const row = [];
    tr.querySelectorAll("input").forEach((inp) => row.push(inp.value.trim() || "0"));
    out.push(row);
  });
  return out;
}

function readB() {
  return readGrid(gridb).map((r) => [r[0]]);
}

function sizeValidate(el) {
  let v = parseInt(el.value, 10);
  if (isNaN(v) || v < 1) v = 1;
  if (v > 6) v = 6;
  el.value = v;
  return v;
}

function currentOpMeta() {
  return OPS[opSel.value] || {};
}

function updateDims() {
  $("dimA").textContent = `(${rowsA.value}×${colsA.value})`;
  $("dimB").textContent = `(${rowsB.value}×${colsB.value})`;
}

function updateExpandOptions() {
  const sel = $("expand");
  const rA = sizeValidate(rowsA);
  const cA = sizeValidate(colsA);
  sel.innerHTML = "";
  for (let i = 1; i <= rA; i++) {
    const o = document.createElement("option");
    o.value = "r" + i;
    o.textContent = `Строка ${i}`;
    if (i === 1) o.selected = true;
    sel.appendChild(o);
  }
  for (let j = 1; j <= cA; j++) {
    const o = document.createElement("option");
    o.value = "c" + j;
    o.textContent = `Столбец ${j}`;
    sel.appendChild(o);
  }
}

function updateSections() {
  const meta = currentOpMeta();
  $("sec-B").classList.toggle("hidden", !meta.b);
  $("sec-k").classList.toggle("hidden", !meta.k);
  $("sec-expand").classList.toggle("hidden", !meta.expand);
  $("sec-b").classList.toggle("hidden", !meta.rhs);
  updateDims();
  if (meta.expand) updateExpandOptions();
}

function rebuild() {
  const rA = sizeValidate(rowsA), cA = sizeValidate(colsA);
  const rB = sizeValidate(rowsB), cB = sizeValidate(colsB);
  buildGrid(gridA, rA, cA, currentA());
  buildGrid(gridB, rB, cB, currentB());
  buildGrid(gridb, rA, 1, currentBvec());
  updateSections();
}

let cache = { A: null, B: null, b: null };

function rebuild() {
  const rA = sizeValidate(rowsA), cA = sizeValidate(colsA);
  const rB = sizeValidate(rowsB), cB = sizeValidate(colsB);
  buildGrid(gridA, rA, cA, cache.A);
  buildGrid(gridB, rB, cB, cache.B);
  buildGrid(gridb, rA, 1, cache.b);
  updateSections();
}

function applyExample() {
  rowsA.value = 3; colsA.value = 3; rowsB.value = 3; colsB.value = 3;
  cache.A = EXAMPLE_A;
  cache.B = EXAMPLE_B;
  cache.b = EXAMPLE_b;
  rebuild();
}

function div(cls) {
  const d = document.createElement("div");
  d.className = cls;
  return d;
}

function buildMatrixBlock(s) {
  const wrap = div("matrix-block" + (s.ans ? " ans" : ""));
  if (s.title) {
    const cap = div("mtitle");
    if (s.title.indexOf("\\") >= 0 && window.katex) {
      katex.render(s.title, cap, { displayMode: false, throwOnError: false });
    } else {
      cap.textContent = s.title;
    }
    wrap.appendChild(cap);
  }
  if (s.caption) {
    const mcap = div("mcap");
    mcap.textContent = s.caption;
    wrap.appendChild(mcap);
  }
  const holder = div("matrix-latex");
  const ltx = matLatex(s.data);
  if (window.katex) katex.render(ltx, holder, { displayMode: false, throwOnError: false });
  else holder.textContent = ltx;
  wrap.appendChild(holder);
  return wrap;
}

function matLatex(data) {
  let spec = null, brace = "pmatrix";
  const rows = data.map((row) => {
    const numeric = row.filter((c) => !(c && c.sep));
    return numeric.join(" & ");
  });
  const first = data[0] || [];
  let sepAt = -1;
  for (const row of data) {
    let cnt = 0;
    for (const c of row) {
      if (c && c.sep) { sepAt = cnt; break; }
      cnt++;
    }
    if (sepAt >= 0) break;
  }
  if (sepAt >= 0) {
    const total = first.filter((c) => !(c && c.sep)).length;
    spec = "c".repeat(sepAt) + "|" + "c".repeat(Math.max(total - sepAt, 0));
    brace = "array";
    return "\\left[\\begin{array}{" + spec + "}" + rows.join(" \\\\ ") + "\\end{array}\\right]";
  }
  return "\\begin{pmatrix}" + rows.join(" \\\\ ") + "\\end{pmatrix}";
}

function renderStepsRaw(box, steps) {
  box.innerHTML = "";
  let stepNum = 0;
  for (const s of steps) {
    if (s.t === "h") {
      const el = div("step-h " + (s.l === 2 ? "l2" : "l1"));
      if (s.l !== 2) {
        stepNum++;
        const n = document.createElement("span");
        n.className = "num";
        n.textContent = "Шаг " + stepNum + ".";
        el.appendChild(n);
      }
      el.appendChild(document.createTextNode(s.s));
      box.appendChild(el);
    } else if (s.t === "p") {
      const el = div("step-p" + (s.ans ? " ans" : ""));
      el.textContent = s.s;
      box.appendChild(el);
    } else if (s.t === "l") {
      const el = div("step-l" + (s.ans ? " ans" : ""));
      if (window.katex) katex.render(s.s, el, { displayMode: false, throwOnError: false });
      else el.textContent = s.s;
      box.appendChild(el);
    } else if (s.t === "m") {
      box.appendChild(buildMatrixBlock(s));
    }
  }
}

function renderSteps(steps) {
  renderStepsRaw($("solution"), steps);
}

function copyText(text, feed) {
  const done = () => {
    if (feed) {
      const prev = feed.textContent;
      feed.textContent = "Скопировано ✔";
      feed.disabled = true;
      setTimeout(() => {
        feed.textContent = prev;
        feed.disabled = false;
      }, 1500);
    }
  };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done).catch(() => legacyCopy(text, done));
  } else {
    legacyCopy(text, done);
  }
}

function legacyCopy(text, done) {
  const ta = document.createElement("textarea");
  ta.value = text;
  ta.style.position = "fixed";
  ta.style.opacity = "0";
  document.body.appendChild(ta);
  ta.select();
  try {
    document.execCommand("copy");
  } catch (e) { /* ignore */ }
  document.body.removeChild(ta);
  done();
}

function renderLatex(res) {
  const box = $("solution");
  box.innerHTML = "";

  const note = div("latex-page-note");
  note.textContent = "Это полный LaTeX-документ. Скопируйте код и скомпилируйте его (например, в Overleaf) либо распечатайте просмотр ниже.";
  box.appendChild(note);

  const bar = div("latex-bar");
  const copyBtn = document.createElement("button");
  copyBtn.textContent = "Скопировать код";
  copyBtn.addEventListener("click", () => copyText(res.latex, copyBtn));
  const printBtn = document.createElement("button");
  printBtn.textContent = "Печать";
  printBtn.addEventListener("click", () => window.print());
  bar.appendChild(copyBtn);
  bar.appendChild(printBtn);
  box.appendChild(bar);

  const codeEl = document.createElement("pre");
  codeEl.className = "latex-code";
  codeEl.textContent = res.latex;
  box.appendChild(codeEl);

  const sec = document.createElement("h2");
  sec.className = "latex-preview-title";
  sec.textContent = "Просмотр решения";
  box.appendChild(sec);

  const prev = div("latex-preview");
  box.appendChild(prev);
  renderPure(prev, res.steps);
}

function renderPure(box, steps) {
  box.innerHTML = "";
  for (const s of steps) {
    if (s.t === "h" || s.t === "p") continue;
    if (s.t === "l") {
      const el = div("step-l");
      const ltx = s.s.replace(/\\text\{[^{}]*\}/g, "").trim();
      if (window.katex) katex.render(ltx, el, { displayMode: false, throwOnError: false });
      else el.textContent = ltx;
      box.appendChild(el);
    } else if (s.t === "m") {
      let ltx = matLatex(s.data || []);
      if (s.title && /[\\_{}^]/.test(s.title)) ltx = s.title + " = " + ltx;
      const el = div("matrix-latex");
      if (window.katex) katex.render(ltx, el, { displayMode: false, throwOnError: false });
      else el.textContent = ltx;
      const wrap = div("matrix-line");
      wrap.appendChild(el);
      box.appendChild(wrap);
    }
  }
}

function renderError(msg) {
  const box = $("solution");
  const err = div("error-box");
  err.textContent = "Ошибка: " + msg;
  box.innerHTML = "";
  box.appendChild(err);
}

let currentMode = "normal";

function setMode(mode) {
  currentMode = mode;
  $("mode-normal").classList.toggle("active", mode === "normal");
  $("mode-latex").classList.toggle("active", mode === "latex");
  $("mode-note").classList.toggle("hidden", mode !== "latex");
}

function solve() {
  const meta = currentOpMeta();
  const A = readGrid(gridA);
  const payload = { op: opSel.value, A: A, mode: currentMode };

  if (meta.b) payload.B = readGrid(gridB);
  if (meta.k) payload.k = $("k").value.trim() || "1";
  if (meta.expand) {
    payload.expand = { kind: $("expand").value[0] === "c" ? "col" : "row", idx: parseInt($("expand").value.slice(1), 10) - 1 };
  }
  if (meta.rhs) {
    payload.b = readGrid(gridb).map((r) => [r[0]]);
  }

  fetch("/api/calc", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  })
    .then((r) => r.json())
    .then((res) => {
      if (res.error) renderError(res.error);
      else if (currentMode === "latex" && res.latex) renderLatex(res);
      else renderSteps(res.steps);
    })
    .catch((err) => renderError(String(err)));
}

opSel.addEventListener("change", () => rebuild());
rowsA.addEventListener("change", () => rebuild());
colsA.addEventListener("change", () => rebuild());
rowsB.addEventListener("change", () => rebuild());
colsB.addEventListener("change", () => rebuild());
$("solve").addEventListener("click", solve);
$("example").addEventListener("click", applyExample);
$("mode-normal").addEventListener("click", () => setMode("normal"));
$("mode-latex").addEventListener("click", () => setMode("latex"));

applyExample();
setMode("normal");
updateSections();