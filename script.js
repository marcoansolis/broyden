const $ = s => document.querySelector(s);
const FN = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'sqrt', 'abs', 'pow', 'min', 'max'];
const COL = ['#0f766e', '#b45309', '#2563eb', '#be185d', '#7c3aed', '#65a30d', '#dc2626', '#0891b2'];
let S = null, timer = null, deb = null;

const esc = s => String(s).replace(/"/g, '&quot;');

function build(f, x) {
    const n = Math.min(8, Math.max(1, +$('#n').value || 1));
    $('#n').value = n;
    $('#eqs').innerHTML = Array.from({ length: n }, (_, i) => `<label>F${i + 1}(x) = 0<input class="eq" value="${esc((f && f[i]) || '')}"></label>`).join('');
    $('#x0s').innerHTML = Array.from({ length: n }, (_, i) => `<label>Valor inicial ${vn(i, n)}<input class="x0" type="number" step="any" value="${x && x[i] !== undefined ? x[i] : 1}"></label>`).join('');
}

const cur = () => ({
    f: [...document.querySelectorAll('.eq')].map(e => e.value),
    x: [...document.querySelectorAll('.x0')].map(e => e.value)
});

function compile(src, n) {
    let s = src.trim();
    if (!s) throw Error('Falta una ecuación.');
    if (s.includes('=')) {
        const p = s.split('=');
        if (p.length !== 2) throw Error('Usa una sola igualdad.');
        s = '(' + p[0] + ')-(' + p[1] + ')';
    }
    if (n <= 3) s = s.replace(/\bx\b/g, 'x1').replace(/\by\b/g, 'x2').replace(/\bz\b/g, 'x3');
    s = s.replace(/\^/g, '**').replace(/\bln\b/g, 'log').replace(/\bpi\b/g, 'PI');
    const vars = Array.from({ length: n }, (_, i) => 'x' + (i + 1));
    let g;
    try {
        g = new Function(...vars, ...FN, 'PI', 'E', 'return (' + s + ');');
    } catch (e) {
        throw Error('Expresión no válida');
    }
    const m = FN.map(k => Math[k]);
    return x => g(...x, ...m, Math.PI, Math.E);
}

const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const nInf = v => Math.max(...v.map(Math.abs));

function solve(A, b) {
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let c = 0; c < n; c++) {
        let p = c;
        for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
        if (Math.abs(M[p][c]) < 1e-14) return null;
        [M[c], M[p]] = [M[p], M[c]];
        for (let r = c + 1; r < n; r++) {
            const m = M[r][c] / M[c][c];
            for (let k = c; k <= n; k++) M[r][k] -= m * M[c][k];
        }
    }
    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let s = M[i][n];
        for (let k = i + 1; k < n; k++) s -= M[i][k] * x[k];
        x[i] = s / M[i][i];
    }
    return x;
}

function jac(F, x) {
    const f0 = F(x), n = x.length, J = f0.map(() => Array(n).fill(0));
    for (let j = 0; j < n; j++) {
        const h = 1e-6 * Math.max(1, Math.abs(x[j])), xp = [...x];
        xp[j] += h;
        const f1 = F(xp);
        for (let i = 0; i < n; i++) J[i][j] = (f1[i] - f0[i]) / h;
    }
    return J;
}

function calc() {
    stop();
    try {
        if ([...document.querySelectorAll('.eq')].some(e => !e.value.trim())) {
            $('#msg').textContent = '';
            return;
        }
        const n = +$('#n').value;
        const fs = [...document.querySelectorAll('.eq')].map(e => compile(e.value, n));
        const F = x => fs.map(f => f(x));
        let x = [...document.querySelectorAll('.x0')].map(e => parseFloat(e.value));
        const tol = parseFloat($('#tol').value.replace(',', '.'));
        const max = +$('#max').value;
        let Fx = F(x);
        let B = $('#jm').value === 'id' ? x.map((_, i) => x.map((_, j) => +(i === j))) : jac(F, x);
        const it = [];
        let pe = null, pp = null;

        for (let k = 0; k <= max; k++) {
            const nF = nInf(Fx);
            const r = { k, x, Fx, nF, B, J: jac(F, x), e: pe, ep: pp };
            it.push(r);
            if (nF < tol || (pe !== null && pe < tol)) break;
            if (k === max) break;
            const s = solve(B, Fx.map(v => -v));
            if (!s) break;
            const xn = x.map((v, i) => v + s[i]), Fn = F(xn);
            const y = Fn.map((v, i) => v - Fx[i]);
            const Bs = B.map(w => dot(w, s));
            const ss = dot(s, s);
            const Bn = B.map((w, i) => w.map((v, j) => v + (y[i] - Bs[i]) * s[j] / ss));
            pe = nInf(s);
            pp = 100 * Math.max(...s.map((v, i) => xn[i] ? Math.abs(v / xn[i]) : Math.abs(v)));
            Object.assign(r, { s, xn, Fn, y, Bn, ss });
            x = xn; Fx = Fn; B = Bn;
        }
        S = { it, fs, n, status: "Proceso completado", ok: true, k: 0 };
        $('#msg').textContent = '';
        setup();
    } catch (e) {
        $('#msg').textContent = e.message;
    }
}

const f = v => {
    if (!isFinite(v)) return String(v);
    const s = v.toFixed(+$('#dec').value);
    return +s === 0 ? s.replace('-', '') : s;
};

const vn = (i, n) => n <= 3 ? 'xyz'[i] : 'x' + (i + 1);
const vec = (v, p) => v.map((a, i) => (p ? (p === 'x' ? vn(i, v.length) : p + (i + 1)) + ' = ' : '') + f(a)).join('<br>');
const mat = M => '<div class="sc"><table class="m">' + M.map(r => '<tr>' + r.map(v => `<td>${f(v)}</td>`).join('') + '</tr>').join('') + '</table></div>';

function plot(items, xr, yr, xl, yl) {
    const W = 520, H = 320, m = { l: 54, r: 12, t: 10, b: 32 };
    const X = v => m.l + (v - xr[0]) / (xr[1] - xr[0]) * (W - m.l - m.r);
    const Y = v => H - m.b - (v - yr[0]) / (yr[1] - yr[0]) * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}">`;
    for (const t of items) {
        if (t.segs) s += `<path d="${t.segs.map(a => `M${X(a[0])} ${Y(a[1])}L${X(a[2])} ${Y(a[3])}`).join('')}" fill="none" stroke="${t.c}" stroke-width="1.6"/>`;
        if (t.line) s += `<polyline fill="none" stroke="${t.c}" stroke-width="1.8"${t.d ? ' stroke-dasharray="6 4"' : ''} points="${t.line.map(p => X(p[0]) + ',' + Y(p[1])).join(' ')}"/>`;
        if (t.dots) t.dots.forEach(p => s += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="${t.r || 3.5}" fill="${t.c}"/>`);
    }
    return s + `</svg>`;
}

function contour(g, xr, yr, N = 70) {
    const xs = i => xr[0] + (xr[1] - xr[0]) * i / N;
    const ys = j => yr[0] + (yr[1] - yr[0]) * j / N;
    const V = [], seg = [];
    for (let i = 0; i <= N; i++) {
        V[i] = [];
        for (let j = 0; j <= N; j++) V[i][j] = g(xs(i), ys(j));
    }
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            const c = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], p = [];
            for (let e = 0; e < 4; e++) {
                const a = c[e], b = c[(e + 1) % 4], va = V[a[0]][a[1]], vb = V[b[0]][b[1]];
                if (va * vb < 0) {
                    const t = va / (va - vb);
                    p.push([xs(a[0]) + (xs(b[0]) - xs(a[0])) * t, ys(a[1]) + (ys(b[1]) - ys(a[1])) * t]);
                }
            }
            if (p.length >= 2) seg.push([...p[0], ...p[1]]);
        }
    }
    return seg;
}

const rng = (v, k = .6) => [
    Math.min(...v) - Math.max(1, (Math.max(...v) - Math.min(...v)) * k),
    Math.max(...v) + Math.max(1, (Math.max(...v) - Math.min(...v)) * k)
];

function setup() {
    const { it, n, fs } = S, K = it.length - 1;
    $('#out').hidden = false;
    $('#sl').max = K;
    if (n === 1) {
        const R = rng(it.map(r => r.x[0])), c = [];
        for (let i = 0; i <= 200; i++) {
            const a = R[0] + (R[1] - R[0]) * i / 200, v = fs[0]([a]);
            if (isFinite(v)) c.push([a, v]);
        }
        S.R = R; S.c = c; S.Y = [-5, 5];
    } else if (n === 2) {
        S.R1 = rng(it.map(r => r.x[0]));
        S.R2 = rng(it.map(r => r.x[1]));
        S.c1 = contour((a, b) => fs[0]([a, b]), S.R1, S.R2);
        S.c2 = contour((a, b) => fs[1]([a, b]), S.R1, S.R2);
    }
    $('#tb').innerHTML = '<tr><th>k</th><th>x⁽ᵏ⁾</th><th>F(x⁽ᵏ⁾)</th><th>‖F‖∞</th><th>ε %</th></tr>' +
        it.map(r => `<tr data-k="${r.k}"><td>${r.k}</td><td>${vec(r.x, 'x')}</td><td>${vec(r.Fx, 'F')}</td><td>${f(r.nF)}</td><td>${r.ep === null ? '—' : f(r.ep)}</td></tr>`).join('');
    view();
}

function view() {
    const { it, n, k } = S, r = it[k], K = it.length - 1, nx = r.xn;
    $('#sl').value = k;
    $('#kl').textContent = `Iteración ${k} de ${K}`;
    if (n === 1) {
        const b = r.B[0][0], [a, c] = S.R, x = r.x[0], y = r.Fx[0], ln = m => [[a, y + m * (a - x)], [c, y + m * (c - x)]];
        $('#ch').innerHTML = plot([{ line: S.c, c: COL[2] }, { line: ln(b), c: COL[0] }, { dots: it.slice(0, k + 1).map(q => [q.x[0], q.Fx[0]]), c: COL[0] }].concat(nx ? [{ dots: [[nx[0], 0]], c: COL[3], r: 5 }] : []), S.R, S.Y, 'x', 'F(x)');
    } else if (n === 2) {
        const P = it.slice(0, k + 1).map(q => [q.x[0], q.x[1]]);
        $('#ch').innerHTML = plot([{ segs: S.c1, c: COL[2] }, { segs: S.c2, c: COL[3] }, { line: P, c: COL[0] }, { dots: P, c: COL[0] }].concat(nx ? [{ dots: [[nx[0], nx[1]]], c: COL[1], r: 5 }] : []), S.R1, S.R2, 'x', 'y');
    }
    $('#p1').innerHTML = vec(r.x, 'x') + '<br><br>' + vec(r.Fx, 'F') + `<br><br>‖F‖∞ = ${f(r.nF)}`;
    $('#p2').innerHTML = `<p class="p"><b>B</b> (Broyden)</p>${mat(r.B)}<p class="p"><b>J</b> (Jacobiana Real)</p>${mat(r.J)}`;
    $('#p3').innerHTML = r.s ? `<p class="p">① Paso s</p><div class="v">${vec(r.s, 's')}</div><p class="p">② Nuevo x</p><div class="v">${vec(r.xn, 'x')}</div><p class="p">③ Actualización B</p>${mat(r.Bn)}` : '';
    document.querySelectorAll('#tb tr[data-k]').forEach(t => t.classList.toggle('cur', +t.dataset.k === k));
}

function go(k) {
    if (!S) return;
    S.k = Math.min(S.it.length - 1, Math.max(0, k));
    view();
}

function stop() {
    clearInterval(timer);
    timer = null;
    $('#pl').textContent = 'Reproducir';
}

$('#pl').onclick = () => {
    if (!S) return;
    if (timer) { stop(); return; }
    $('#pl').textContent = 'Pausar';
    timer = setInterval(() => {
        if (S.k >= S.it.length - 1) { stop(); return; }
        go(S.k + 1);
    }, 900);
};

$('#pv').onclick = () => { stop(); go(S.k - 1); };
$('#nx').onclick = () => { stop(); go(S.k + 1); };
$('#sl').oninput = e => { stop(); go(+e.target.value); };
$('#tb').onclick = e => { const t = e.target.closest('tr[data-k]'); if (t) { stop(); go(+t.dataset.k); } };
$('#n').onchange = () => { const c = cur(); build(c.f, c.x); calc(); };
$('#setup').addEventListener('input', e => { if (e.target.id === 'n') return; clearTimeout(deb); deb = setTimeout(calc, 300); });

build(['x^2+y^2-4', 'x^2-y-1'], [1, 1]);
calc();
