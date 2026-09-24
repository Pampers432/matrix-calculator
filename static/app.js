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
    cap.textContent = s.title;
    wrap.appendChild(cap);
  }
  if (s.caption) {
    const mcap = div("mcap");
    mcap.textContent = s.caption;
    wrap.appendChild(mcap);
  }
  const tb = div("paren matrixbox");
  const table = document.createElement("table");
  table.className = "matrix";
  for (const row of s.data) {
    const tr = document.createElement("tr");
    for (const c of row) {
      const td = document.createElement("td");
      if (c && typeof c === "object" && c.sep) {
        td.className = "sep";
        td.textContent = "|";
      } else {
        td.innerHTML = c;
      }
      tr.appendChild(td);
    }
    table.appendChild(tr);
  }
  tb.appendChild(table);
  wrap.appendChild(tb);
  return wrap;
}

function renderSteps(steps) {
  const box = $("solution");
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
      el.innerHTML = s.s;
      box.appendChild(el);
    } else if (s.t === "m") {
      box.appendChild(buildMatrixBlock(s));
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

function solve() {
  const meta = currentOpMeta();
  const A = readGrid(gridA);
  const payload = { op: opSel.value, A: A };

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

applyExample();
updateSections();