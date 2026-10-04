"""Revision 3, topic `numbers`: traceability audit of every number in the revision-2 manuscript and supplement.

Read-only analysis (no manuscript source, script or CSV is modified). Usage (from the repository root or analysis/):

    python3 analysis/rev3_numbers_audit.py [--rev2-out <dir with regenerated rev2_summary.json / rev2_tables.tex>]
                                           [--json analysis/rev3_numbers_audit.json]

If --rev2-out does not contain rev2_summary.json, `analysis/rev2_analysis.py --out-dir <dir>` is run first (about 15 s).

Parts (numbering as in the task):
 1. CHECK comments: every "% CHECK..." comment of SWEVO_rev2/latex_source/optA/swevo_front.tex, optA/sw/*.tex
    (incl. supp_theory.tex) and SWEVO_supplement.tex; each \\N... macro the comment names (with or without the
    backslash) is resolved to its value in the REPOSITORY's generated macro files (analysis/mpce_numbers*.tex) and
    searched, as printed, in the paragraph the comment belongs to (lines before the comment up to the previous blank
    line or comment block; an inline comment includes its own line).
 2. Literal numbers: every number in the non-comment text of the revision-2 main text (front matter + sections,
    including inline tables) and of the supplement prose is traced
      (a) positionally: the revision-2 files are aligned token by token with the repository's macro versions of the
          same files (optA/swevo_front.tex, optA/sw/*.tex, SWEVO_supplement.tex), macros expanded with the
          repository values; a number at the position of a macro is "anchored" (equal or different);
      (b) by value: macro values, generated tables (analysis/mpce_tab_*.tex, mpce_supp_*.tex, mpce_supplementary.tex),
          regenerated rev2_summary.json / rev2_tables.tex, SWEVO_rev2/validation/*.json|*.tex|*.csv,
          SWEVO_rev2/latex_source/analysis/archive_reliability_tables.tex, the repository summaries
          analysis/mpce_summary*.json (primary pool); other analysis/*.tex|*.json and small non-run CSVs (secondary).
    Classes: anchored_equal, anchored_different, traced (primary pool), traced_secondary, structural (design
    constants, small integers, years ...), untraced.
 3. Hand-edited generated files: SWEVO_rev2/latex_source/analysis/*.tex vs analysis/*.tex, every changed number;
    the per-evaluation times are recomputed from the run records (Seconds / Calls of the 6,030-evaluation runs of the
    68 cases, random initialization), median (definition in mpce_results.py) and mean (Table S-cost).
Check-suite results (part 4) are collected by the lead/agent from the logs named in the report.
"""
import os, re, sys, json, glob, math, argparse, subprocess, difflib, csv
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
R2 = os.path.join(ROOT, "SWEVO_rev2", "latex_source")
VAL = os.path.join(ROOT, "SWEVO_rev2", "validation")

MAIN_FILES = ["optA/swevo_front.tex"] + [f"optA/sw/{b}.tex" for b in (
    "02_intro", "02b_related", "03_model", "04_methods", "05_setup", "06_results", "07_ablation", "08_beyond",
    "09_robust", "10_limits_concl")]
SUPP_FILES = ["SWEVO_supplement.tex", "optA/sw/supp_theory.tex"]
# repository counterparts (macro versions) of the revision-2 files, for the positional alignment
REPO_PAIR = {f: f for f in MAIN_FILES + SUPP_FILES}

# ------------------------------------------------------------------------------------------------ macros
NEWCMD = re.compile(r"\\(?:newcommand|renewcommand|providecommand)\*?\{\\(N[A-Za-z]+)\}\{(.*)\}\s*$")


def load_macros(files):
    d = {}
    for fn in files:
        for line in open(fn, encoding="utf-8"):
            m = NEWCMD.match(strip_comment(line).strip())
            if m and m.group(1) not in d:
                d[m.group(1)] = m.group(2)
    return d


MACRO_USE = re.compile(r"\\(N[A-Z][A-Za-z]*)(?:\{\})?")


def expand(s, M, depth=0):
    if depth > 6:
        return s
    out = MACRO_USE.sub(lambda m: M.get(m.group(1), m.group(0)), s)
    return out if out == s else expand(out, M, depth + 1)


# ------------------------------------------------------------------------------------------------ LaTeX -> plain
def strip_comment(line):
    m = re.search(r"(?<!\\)%", line)
    return line if m is None else line[:m.start()]


def remove_cmd_args(s, cmds):
    """Remove \\cmd[opt]{arg}{arg}... (balanced braces) for the commands given (name -> number of {} args)."""
    out, i = [], 0
    pat = re.compile(r"\\(%s)\b\*?" % "|".join(sorted(cmds, key=len, reverse=True)))
    while True:
        m = pat.search(s, i)
        if not m:
            out.append(s[i:]); break
        out.append(s[i:m.start()])
        j = m.end(); nargs = cmds[m.group(1)]
        while j < len(s) and s[j] in " \t":
            j += 1
        if j < len(s) and s[j] == "[":
            k = s.find("]", j)
            j = k + 1 if k > 0 else j
        for _ in range(nargs):
            while j < len(s) and s[j] in " \t\n":
                j += 1
            if j < len(s) and s[j] == "{":
                depth = 0
                for k in range(j, len(s)):
                    depth += (s[k] == "{") - (s[k] == "}")
                    if depth == 0:
                        break
                j = k + 1
        out.append(" ")
        i = j
    return "".join(out)


DROP_CMDS = {"label": 1, "ref": 1, "eqref": 1, "cite": 1, "citep": 1, "citet": 1, "input": 1, "includegraphics": 1,
             "pipefig": 3, "SuppTable": 1, "SuppTableOpt": 1, "SuppOmit": 1, "url": 1, "texttt": 1, "hspace": 1,
             "vspace": 1, "setlength": 2, "rule": 2, "cmidrule": 1, "cline": 1, "begin": 1, "end": 1,
             "IfFileExists": 1, "CaptureSuppTables": 1, "IfSuppTable": 1, "tabcolsep": 0, "arraystretch": 0,
             "resizebox": 2, "addlinespace": 0, "nolinenumbers": 0, "href": 1, "bibitem": 1, "newblock": 0,
             "doi": 1, "fbox": 0, "parbox": 1, "multirow": 2}


def tex_to_plain(s):
    s = remove_cmd_args(s, DROP_CMDS)
    s = re.sub(r"\\multicolumn\{\d+\}\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", " ", s)
    s = re.sub(r"(?:p|m|b)\{[\d.]*\\(?:textwidth|columnwidth|linewidth)\}", " ", s)
    s = re.sub(r"[\d.]+\\(?:textwidth|columnwidth|linewidth|baselineskip)", " ", s)
    s = re.sub(r"\{(?:@\{\}|[lcrX|]|>\{[^}]*\}|p\{[^}]*\})+\}", " ", s)            # tabular column specs
    s = re.sub(r"\\\\\[[^\]]*\]", " ", s)
    s = re.sub(r"\b[\d.]+\s*(?:pt|cm|mm|em|ex|bp)\b", " ", s)
    s = s.replace("\\ensuremath", "")
    s = re.sub(r"\{,\}", "", s)
    s = re.sub(r"(\d)\\,(\d{3})", r"\1\2", s)
    s = re.sub(r"\$?\^\{?\\circ\}?\$?|\\textdegree(?:\{\})?", "°", s)
    s = s.replace("\\%", "%")
    s = re.sub(r"(\d+(?:\.\d+)?)\s*(?:\\times|\\cdot)\s*\{?\s*10\s*\}?\^\{\s*(?:\{-\}|-|\\textminus\s*)\s*(\d+)\s*\}",
               r"\1e-\2", s)
    s = re.sub(r"(\d+(?:\.\d+)?)\s*(?:\\times|\\cdot)\s*\{?\s*10\s*\}?\^\{?\s*\+?(\d+)\s*\}?", r"\1e\2", s)
    s = re.sub(r"(?<![\d.])10\^\{\s*(?:\{-\}|-)\s*(\d+)\s*\}", r"1e-\1", s)
    s = re.sub(r"(?<![\d.])10\^\{?(\d+)\}?", r"1e\1", s)
    s = s.replace("--", " – ")
    s = re.sub(r"\$\s*-\s*\$|\{-\}|\\textminus(?:\{\})?|−|\\mbox\{-\}", "-", s)
    s = re.sub(r"\$\s*\+\s*\$|\{\+\}", "+", s)
    s = re.sub(r"\\pm", " ± ", s)
    s = re.sub(r"[_^]\{[^{}]*\}", " ", s)          # sub/superscripts (indices, not values)
    s = re.sub(r"[_^]\\?[A-Za-z0-9]", " ", s)
    s = re.sub(r"\\[A-Za-z]+\*?", " ", s)          # remaining commands
    s = s.replace("~", " ").replace("$", " ").replace("{", " ").replace("}", " ").replace("&", " & ")
    s = re.sub(r"(\d),(\d{3})(?![\d])", r"\1\2", s)
    s = re.sub(r"(\d),(\d{3})(?![\d])", r"\1\2", s)
    return s


NUM = re.compile(r"(?P<sign>[+\-±]?)(?P<num>\d+(?:\.\d+)?)(?:e(?P<exp>[+-]?\d+))?")


def extract_numbers(plain):
    """[(start, signed string, unsigned key, value, decimals, sig digits, pct)] of the numbers in plain text."""
    res = []
    for m in NUM.finditer(plain):
        st = m.start("num")
        prev = plain[st - 1] if st > 0 else " "
        prev_sign = plain[m.start() - 1] if m.start() > 0 else " "
        nxt = plain[m.end():m.end() + 1]
        if prev.isalpha() or prev == "." or (prev.isdigit() and not m.group("sign")):
            continue                                   # part of a word (IEA37, S1, V80) or of another number
        if nxt.isalpha() and nxt not in "e":           # 4D, 5D, 1st, 12pt
            continue
        if nxt.isalpha() and nxt == "e" and not m.group("exp"):
            continue
        sign = m.group("sign")
        if sign and (prev_sign.isalnum() or prev_sign in ")]}"):
            sign = ""                                  # hyphen / range, not a sign
        num, exp = m.group("num"), m.group("exp")
        dec = len(num.split(".")[1]) if "." in num else 0
        digits = num.replace(".", "").lstrip("0")
        sig = max(len(digits), 1)
        key = num if exp is None else f"{num}e{int(exp)}"
        try:
            val = float(num) * (10 ** int(exp) if exp else 1)
        except OverflowError:
            val = float("inf")
        if sign == "-":
            val = -val
        pct = plain[m.end():m.end() + 2].lstrip().startswith("%")
        res.append(dict(pos=m.start(), sign=sign if sign in "+-" else "", key=key, value=val, dec=dec, sig=sig,
                        pct=pct, exp=exp is not None, end=m.end()))
    return res


# ------------------------------------------------------------------------------------------------ source pools
def walk_json(x, p=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from walk_json(v, f"{p}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from walk_json(v, f"{p}/{i}")
    else:
        yield p, x


class Pool:
    """Value index: unsigned printed key -> [source descriptions]."""

    def __init__(self):
        self.idx = defaultdict(list)

    def add_key(self, key, src):
        lst = self.idx[key]
        if len(lst) < 6 and src not in lst:
            lst.append(src)

    def add_text(self, text, src_prefix, per_line=True):
        for no, line in enumerate(text.splitlines(), 1):
            line = strip_comment(line)
            for n in extract_numbers(tex_to_plain(line)):
                self.add_key(n["key"], f"{src_prefix}:{no}")

    def add_float(self, x, src):
        if isinstance(x, bool) or x is None:
            return
        if isinstance(x, str):
            try:
                x = float(x)
            except ValueError:
                for n in extract_numbers(x):
                    self.add_key(n["key"], src)
                return
        if not isinstance(x, (int, float)) or math.isnan(x) or math.isinf(x):
            return
        for y, tag in ((abs(x), ""), (abs(x) * 100, "x100")):
            s = src + (f" {tag}" if tag else "")
            if abs(y - round(y)) < 1e-9 and y < 1e9:
                self.add_key(str(int(round(y))), s)
            for k in range(1, 5):
                self.add_key(f"{y:.{k}f}", s)
            if 0 < y < 0.01 or y >= 1e5:
                for k in range(0, 3):
                    m, e = f"{y:.{k}e}".split("e")
                    self.add_key(f"{m}e{int(e)}", s)

    def get(self, key):
        return self.idx.get(key, [])


def build_pools(rev2_out):
    prim, sec, macro, pj = Pool(), Pool(), Pool(), Pool()     # prim: printed (tex), pj: json / csv values
    r2 = Pool()                                               # revision-2 sources only (rev2 outputs, validation, archive)
    repo_macro_files = sorted(glob.glob(os.path.join(HERE, "mpce_numbers*.tex")))
    M = load_macros(repo_macro_files)
    for name, val in M.items():
        for n in extract_numbers(tex_to_plain(expand(val, M))):
            macro.add_key(n["key"], "\\" + name)
    gen = sorted(glob.glob(os.path.join(HERE, "mpce_tab_*.tex")) + glob.glob(os.path.join(HERE, "mpce_supp*.tex")))
    for fn in gen:
        prim.add_text(open(fn, encoding="utf-8").read(), os.path.relpath(fn, ROOT))
    for fn in [os.path.join(rev2_out, "rev2_tables.tex"), os.path.join(R2, "analysis", "archive_reliability_tables.tex")] \
            + sorted(glob.glob(os.path.join(VAL, "*.tex"))):
        if os.path.exists(fn):
            lab = "rev2_tables.tex(regenerated)" if fn.startswith(rev2_out) else os.path.relpath(fn, ROOT)
            prim.add_text(open(fn, encoding="utf-8").read(), lab)
            r2.add_text(open(fn, encoding="utf-8").read(), lab)
    jsons = [(os.path.join(rev2_out, "rev2_summary.json"), "rev2_summary.json(regenerated)")]
    jsons += [(f, os.path.relpath(f, ROOT)) for f in sorted(glob.glob(os.path.join(VAL, "*.json")))
              if not f.endswith("raw_data_checksums.json")]
    jsons += [(f, os.path.relpath(f, ROOT)) for f in sorted(glob.glob(os.path.join(HERE, "mpce_summary*.json")))]
    for fn, lab in jsons:
        for p, x in walk_json(json.load(open(fn))):
            pj.add_float(x, f"{lab}{p}")
            if "mpce_summary" not in lab:
                r2.add_float(x, f"{lab}{p}")
    for fn in sorted(glob.glob(os.path.join(VAL, "*.csv"))):
        for i, row in enumerate(csv.DictReader(open(fn)), 2):
            for k, v in row.items():
                pj.add_float(v, f"{os.path.relpath(fn, ROOT)}:{i}:{k}")
                r2.add_float(v, f"{os.path.relpath(fn, ROOT)}:{i}:{k}")
    # secondary: other generated tex / json of the repository, small non-run csv
    for fn in sorted(glob.glob(os.path.join(HERE, "*.tex"))):
        b = os.path.basename(fn)
        if b.startswith(("mpce_tab_", "mpce_supp", "mpce_numbers", "rev3_")):
            continue
        sec.add_text(open(fn, encoding="utf-8", errors="replace").read(), os.path.relpath(fn, ROOT))
    for fn in sorted(glob.glob(os.path.join(HERE, "*.json"))):
        b = os.path.basename(fn)
        if b.startswith(("mpce_summary", "rev3_")):
            continue
        try:
            for p, x in walk_json(json.load(open(fn))):
                sec.add_float(x, f"{os.path.relpath(fn, ROOT)}{p}")
        except Exception:
            pass
    for fn in sorted(glob.glob(os.path.join(HERE, "*.csv"))):
        b = os.path.basename(fn)
        if re.search(r"_s\d+of\d+\.csv$", b) or os.path.getsize(fn) > 400_000 or b.startswith("rev3_"):
            continue
        try:
            for i, row in enumerate(csv.DictReader(open(fn, encoding="utf-8", errors="replace")), 2):
                for k, v in row.items():
                    if k and k not in ("Coordinates", "Curve"):
                        sec.add_float(v, f"{os.path.relpath(fn, ROOT)}:{k}")
        except Exception:
            pass
    return M, macro, prim, pj, sec, r2


def density(n, *pools):
    """How many of the 10 neighbours v +- j * 10^-dec (j = 1..5) of a printed number are also in the pools: a value
    found in a pool that also holds most of its neighbours is not specifically traced."""
    if n["exp"]:
        return 0
    v, d = abs(n["value"]), n["dec"]
    c = 0
    for j in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5):
        w = v + j * 10 ** -d
        if w < 0:
            continue
        k = f"{w:.{d}f}" if d else str(int(round(w)))
        c += any(p.get(k) for p in pools)
    return c


# ------------------------------------------------------------------------------------------------ part 1
def check_blocks(fn):
    """[(comment line numbers, comment text, paragraph line numbers)] of the CHECK comments of a file."""
    lines = open(fn, encoding="utf-8").read().split("\n")
    is_com = [bool(re.match(r"\s*%", l)) for l in lines]
    has_check = [("%" in l and "CHECK" in l[l.index("%"):]) if re.search(r"(?<!\\)%", l) else False for l in lines]
    out, i = [], 0
    n = len(lines)
    while i < n:
        if is_com[i]:
            j = i
            while j + 1 < n and is_com[j + 1]:
                j += 1
            if any(has_check[i:j + 1]):
                par = []
                k = i - 1
                while k >= 0 and lines[k].strip() and not is_com[k]:
                    par.insert(0, k); k -= 1
                for c in range(i, j + 1):
                    if has_check[c]:
                        out.append(dict(com_lines=[c + 1], com=lines[c], par=[p + 1 for p in par], inline=False))
            i = j + 1
            continue
        if has_check[i]:                               # inline comment after text
            par = [i]
            k = i - 1
            while k >= 0 and lines[k].strip() and not is_com[k]:
                par.insert(0, k); k -= 1
            out.append(dict(com_lines=[i + 1], com=lines[i][re.search(r"(?<!\\)%", lines[i]).start():],
                            par=[p + 1 for p in par], inline=True))
        i += 1
    return lines, out


def norm_words(s):
    s = tex_to_plain(s)
    s = s.replace("‐", "-")
    return re.sub(r"\s+", " ", s).strip().lower()


def _line_keys(line, M):
    return Counter((n["sign"] if n["sign"] == "-" else "") + n["key"]
                   for n in extract_numbers(tex_to_plain(expand(strip_comment(line), M))))


def locate_missing(lines, b, missing, M_rev2, rel, nm):
    """Where the values missing from the comment's paragraph are printed instead: in a neighbouring paragraph
    (within 12 lines before / 3 after the comment: the comment has become detached from its paragraph), elsewhere in
    the same file, or nowhere (the comment documents a value that the text does not print). Also whether the
    repository (macro) version of the file printed the macro in its text."""
    c0 = b["com_lines"][0] - 1
    near_lines, else_lines = [], []
    for i, l in enumerate(lines):
        if re.match(r"\s*%", l):
            continue
        k = _line_keys(l, M_rev2)
        if all(k.get(m) for m in missing):
            (near_lines if c0 - 12 <= i <= c0 + 3 else else_lines).append(i + 1)
    repo_fn = os.path.join(ROOT, REPO_PAIR.get(rel, rel))
    repo_used = []
    if os.path.exists(repo_fn):
        for i, l in enumerate(open(repo_fn, encoding="utf-8").read().split("\n"), 1):
            if re.search(r"\\%s(?![A-Za-z])" % nm, strip_comment(l)):
                repo_used.append(i)
    where = "adjacent_paragraph" if near_lines else ("elsewhere_in_file" if else_lines else "not_printed_in_file")
    return dict(located=where, located_lines=(near_lines or else_lines)[:6], repo_text_macro_lines=repo_used[:6])


def part1(M_repo, M_rev2):
    res = []
    for rel in MAIN_FILES + SUPP_FILES:
        fn = os.path.join(R2, rel)
        lines, blocks = check_blocks(fn)
        for b in blocks:
            com = b["com"]
            names = []
            for m in re.finditer(r"\\?(N[A-Z][A-Za-z]+)", com):
                nm = m.group(1)
                if nm in M_repo and nm not in names:
                    names.append(nm)
                elif m.group(0).startswith("\\") and nm not in names and nm not in M_repo:
                    names.append(nm)
            if not names:
                continue
            par_text = "\n".join(strip_comment(lines[p - 1]) for p in b["par"])
            # macros still used in the paragraph print the REV2 values (supplement); expand them for the search
            par_exp = expand(par_text, M_rev2)
            par_plain = tex_to_plain(par_exp)
            par_nums = extract_numbers(par_plain)
            keys = Counter((n["sign"] if n["sign"] == "-" else "") + n["key"] for n in par_nums)
            ukeys = Counter(n["key"] for n in par_nums)
            for nm in names:
                r = dict(file=rel, comment_line=b["com_lines"][0], paragraph_lines=[b["par"][0], b["par"][-1]] if b["par"] else [],
                         macro=nm, comment=com.strip()[:300])
                if nm not in M_repo:
                    r.update(status="undefined_macro", repo_value=None); res.append(r); continue
                raw = M_repo[nm]
                val = expand(raw, M_repo)
                r["repo_value"] = raw
                if nm in M_rev2 and M_rev2[nm] != raw:
                    r["rev2_macro_file_value"] = M_rev2[nm]
                r["macro_used_in_paragraph"] = bool(re.search(r"\\%s(?![A-Za-z])" % nm, par_text))
                if not b["par"]:
                    r.update(status="no_paragraph"); res.append(r); continue
                vn = extract_numbers(tex_to_plain(val))
                if vn:
                    missing, signless = [], []
                    for n in vn:
                        k = (n["sign"] if n["sign"] == "-" else "") + n["key"]
                        if keys.get(k):
                            continue
                        if ukeys.get(n["key"]):
                            signless.append(k); continue
                        missing.append(k)
                    if not missing and not signless:
                        r["status"] = "match"
                    elif not missing:
                        r.update(status="match_sign_differs", detail=signless)
                    else:
                        # near values in the paragraph (same magnitude) to show what the text has instead
                        near = []
                        for k in missing:
                            try:
                                v = float(k)
                            except ValueError:
                                continue
                            for n in par_nums:
                                if v != 0 and abs(abs(n["value"]) - abs(v)) <= max(0.25 * abs(v), 10 ** -n["dec"]):
                                    near.append((n["sign"] if n["sign"] == "-" else "") + n["key"])
                        r.update(status="MISMATCH", missing=missing, text_near_values=sorted(set(near))[:8])
                        r.update(locate_missing(lines, b, missing, M_rev2, rel, nm))
                else:
                    w = norm_words(val)
                    r["status"] = "match_text" if w and w in norm_words(par_exp) else "text_value_not_found"
                    r["value_plain"] = w
                res.append(r)
    return res


# ------------------------------------------------------------------------------------------------ part 2a (alignment)
TOK = re.compile(r"\d+(?:\.\d+)?(?:e[+-]?\d+)?|[A-Za-z]+|[^\sA-Za-z\d]")


def tokens_of(plain):
    return [(m.group(0), m.start()) for m in TOK.finditer(plain)]


def file_tokens(text, M, track_macros):
    """Tokens of a whole file (comments stripped, macros expanded with M). If track_macros, every token carries the
    macro name it comes from (None for literal text). Returns tokens, owners, and line number of each token."""
    toks, owners, lines_of = [], [], []
    for no, line in enumerate(text.split("\n"), 1):
        line = strip_comment(line)
        if not line.strip():
            continue
        pieces = []
        pos = 0
        for m in MACRO_USE.finditer(line):
            if m.start() > pos:
                pieces.append((line[pos:m.start()], None))
            nm = m.group(1)
            pieces.append((expand(m.group(0), M), nm if (track_macros and nm in M) else None))
            pos = m.end()
        pieces.append((line[pos:], None))
        if not track_macros:
            pieces = [("".join(p for p, _ in pieces), None)]
        for txt, owner in pieces:
            for t, _ in tokens_of(tex_to_plain(txt)):
                toks.append(t); owners.append(owner); lines_of.append(no)
    return toks, owners, lines_of


def align(rel, M_repo, M_rev2):
    old_fn, new_fn = os.path.join(ROOT, REPO_PAIR[rel]), os.path.join(R2, rel)
    if not os.path.exists(old_fn):
        return None
    ot, oo, ol = file_tokens(open(old_fn, encoding="utf-8").read(), M_repo, True)
    # the revision-2 supplement still uses macros: expand them with the REV2 macro files (what is printed)
    nt, no_, nl = file_tokens(open(new_fn, encoding="utf-8").read(), M_rev2, True)
    sm = difflib.SequenceMatcher(None, ot, nt, autojunk=False)
    ops = sm.get_opcodes()
    new_owner = [None] * len(nt)        # for new tokens: ("equal", macro) if aligned to a macro token
    new_tag = ["new"] * len(nt)         # "repo": token aligned (equal) with the repository text, "new": revision-2 text
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(j1, j2):
                new_tag[k] = "repo"
    # group old macro occurrences: consecutive tokens with the same owner and same line
    occ = []
    i = 0
    while i < len(ot):
        if oo[i]:
            j = i
            while j + 1 < len(ot) and oo[j + 1] == oo[i] and ol[j + 1] == ol[i]:
                j += 1
            occ.append((i, j, oo[i]))
            i = j + 1
        else:
            i += 1
    # map old index -> (tag, new index or block)
    o2n = {}
    for tag, i1, i2, j1, j2 in ops:
        for k in range(i1, i2):
            o2n[k] = (tag, i1, i2, j1, j2)
    results = []
    for (a, b, nm) in occ:
        tags = {o2n[k][0] for k in range(a, b + 1)}
        old_txt = " ".join(ot[a:b + 1])
        if tags == {"equal"}:
            for k in range(a, b + 1):
                tag, i1, i2, j1, j2 = o2n[k]
                new_owner[j1 + (k - i1)] = nm
            results.append(dict(macro=nm, status="equal", old=old_txt, old_line=ol[a],
                                new_line=nl[o2n[a][3] + (a - o2n[a][1])]))
            continue
        if tags == {"delete"}:
            results.append(dict(macro=nm, status="removed", old=old_txt, old_line=ol[a]))
            continue
        # some tokens replaced: show the new block(s)
        blocks = sorted({o2n[k][1:] for k in range(a, b + 1) if o2n[k][0] != "equal"})
        new_txt = " | ".join(" ".join(nt[j1:j2]) for (_, _, j1, j2) in blocks)
        old_blk = " | ".join(" ".join(ot[i1:i2]) for (i1, i2, _, _) in blocks)
        nline = nl[blocks[0][2]] if blocks and blocks[0][2] < len(nl) else None
        val_nums = [t for t in ot[a:b + 1] if re.match(r"\d", t)]
        new_nums = set(t for (_, _, j1, j2) in blocks for t in nt[j1:j2] if re.match(r"\d", t))
        # context: 8 new tokens before / after
        if blocks:
            j1, j2 = blocks[0][2], blocks[-1][3]
            ctx = " ".join(nt[max(0, j1 - 10):j1]) + " [[" + " ".join(nt[j1:j2]) + "]] " + " ".join(nt[j2:j2 + 10])
        else:
            ctx = ""
        status = "replaced_value_kept" if val_nums and all(v in new_nums for v in val_nums) else "REPLACED"
        if not val_nums:
            status = "replaced_text_macro"
        results.append(dict(macro=nm, status=status, old=old_txt, old_block=old_blk[:300], new_block=new_txt[:300],
                            old_line=ol[a], new_line=nline, context=ctx[:400]))
    return results, nt, new_owner, nl, new_tag


# ------------------------------------------------------------------------------------------------ part 2b (all numbers)
YEARS = set(range(1950, 2031))
STRUCT_INT = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 24, 25, 30, 36, 40, 45, 50, 60, 64, 68,
              80, 90, 100, 120, 150, 180, 200, 250, 300, 308, 360, 385, 400, 462, 500, 600, 750, 1000, 1300, 2000,
              3015, 6030, 30030, 120030, 3000, 6000, 30000, 120000, 1000000}
STRUCT_DEC = {"0.5", "0.05", "0.95", "0.9", "0.75", "0.25", "1.5", "2.5", "0.7", "0.7298", "1.49618", "1.496", "2.0",
              "0.01", "0.1", "0.2", "0.04", "0.075", "0.3", "0.6", "0.8", "0.4", "1.0", "0.0", "0.001"}


SKIP_LINE = re.compile(r"\\(?:address|ead|author|cortext|fntext|tnotetext|thanks|bibitem)\b|postcode")


def scan_numbers(rel, M_rev2, pools, align_info):
    M_unused, macro, prim, pj, sec, r2 = pools
    inl = {}                                         # line -> verified inlined generated table
    for t in INLINED:
        if t["file"] == rel:
            for l in range(t["lines"][0], t["lines"][1] + 1):
                inl[l] = t
    fn = os.path.join(R2, rel)
    text = open(fn, encoding="utf-8").read()
    out = []
    nt, new_owner, nl, ntag = align_info[1:5] if align_info else (None, None, None, None)
    anchored = defaultdict(list)    # line -> [(token, macro)]
    origin = defaultdict(lambda: defaultdict(list))   # line -> token -> [origin tags]
    if nt is not None:
        for t, o, l, g in zip(nt, new_owner, nl, ntag):
            if o and re.match(r"\d", t):
                anchored[l].append((t, o))
            if re.match(r"\d", t):
                origin[l][t].append("macro" if o else g)
    in_bib = False
    for no, line in enumerate(text.split("\n"), 1):
        code = strip_comment(line)
        if "\\begin{thebibliography}" in code:
            in_bib = True
        if "\\end{thebibliography}" in code:
            in_bib = False
            continue
        if not code.strip() or in_bib or SKIP_LINE.search(code):
            continue
        plain = tex_to_plain(expand(code, M_rev2))
        used_macros = [m.group(1) for m in MACRO_USE.finditer(code)]
        anc = Counter(t for t, _ in anchored.get(no, []))
        anc_m = defaultdict(list)
        for t, o in anchored.get(no, []):
            anc_m[t].append(o)
        for n in extract_numbers(plain):
            k = n["key"]
            ctx = plain[max(0, n["pos"] - 70):n["end"] + 50].replace("\n", " ")
            ctx = re.sub(r"\s+", " ", ctx).strip()
            rec = dict(file=rel, line=no, text=(n["sign"] if n["sign"] == "-" else "") + k, context=ctx)
            mvals = {u: [x["key"] for x in extract_numbers(tex_to_plain(expand("\\" + u, M_rev2)))] for u in used_macros}
            if no in inl:
                rec.update(cls="inlined_generated_table" if inl[no]["identical_numbers"] else "inlined_table_modified",
                           source=[inl[no]["source"]], specific=True)
                if not inl[no]["identical_numbers"]:
                    rec["in_rev2_sources"] = r2.get(k)[:3]
            elif anc.get(k):
                anc[k] -= 1
                rec.update(cls="anchored_equal", source=["\\" + anc_m[k].pop(0)], specific=True)
            elif any(k in v for v in mvals.values()):
                rec.update(cls="macro_in_rev2_text", source=["\\" + u for u, v in mvals.items() if k in v][:3], specific=True)
            elif macro.get(k) or prim.get(k):
                rec.update(cls="traced", source=(macro.get(k) + prim.get(k))[:4], kind="printed",
                           specific=density(n, macro, prim) <= 2)
            elif pj.get(k):
                rec.update(cls="traced", source=pj.get(k)[:4], kind="json", specific=density(n, pj) <= 2)
            elif sec.get(k):
                rec.update(cls="traced_secondary", source=sec.get(k)[:3], specific=density(n, sec) <= 2)
            else:
                v = abs(n["value"])
                struct = (not n["exp"] and n["dec"] == 0 and (v in STRUCT_INT or int(v) in YEARS or v <= 12)) \
                    or k in STRUCT_DEC
                rec.update(cls="structural" if struct else "untraced", source=[], specific=False)
                if not struct:
                    rec["near"] = near_values(n, macro, prim)
            rec["weak"] = n["sig"] <= 2 and not n["exp"]
            if rec["cls"] in ("traced", "traced_secondary", "untraced") and r2.get(k):
                rec["rev2_source"] = r2.get(k)[:3]
                rec["rev2_specific"] = density(n, r2) <= 2
            og = origin.get(no, {}).get(k) if nt is not None else None
            rec["origin"] = og.pop(0) if og else ("unaligned" if nt is not None else "no_repo_counterpart")
            out.append(rec)
    return out


_NEAR_CACHE = {}


def near_values(n, macro, prim):
    """Values of the primary pools within one unit of the last printed digit (possible rounding / typing slips)."""
    if n["exp"] or n["sig"] < 3:
        return []
    v, d = abs(n["value"]), n["dec"]
    cands = []
    for delta in (-1, 1):
        s = f"{v + delta * 10 ** -d:.{d}f}"
        for src in macro.get(s) + prim.get(s):
            cands.append(f"{s} ({src})")
    return cands[:4]


# ------------------------------------------------------------------------------------------------ part 3
def part3():
    changed = []
    for fn in sorted(glob.glob(os.path.join(R2, "analysis", "*.tex"))):
        b = os.path.basename(fn)
        rp = os.path.join(HERE, b)
        if not os.path.exists(rp):
            changed.append(dict(file=b, status="not in repository (new revision-2 file)"))
            continue
        a_l = open(rp, encoding="utf-8").read().split("\n")
        b_l = open(fn, encoding="utf-8").read().split("\n")
        for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a_l, b_l, autojunk=False).get_opcodes():
            if tag == "equal":
                continue
            for k in range(max(i2 - i1, j2 - j1)):
                ol = a_l[i1 + k] if i1 + k < i2 else ""
                nl_ = b_l[j1 + k] if j1 + k < j2 else ""
                on = [x["key"] for x in extract_numbers(tex_to_plain(ol))]
                nn = [x["key"] for x in extract_numbers(tex_to_plain(nl_))]
                rec = dict(file=b, repo_line=i1 + k + 1, rev2_line=j1 + k + 1,
                           numbers_changed=on != nn, repo_numbers=on if on != nn else [], rev2_numbers=nn if on != nn else [])
                m1, m2 = NEWCMD.match(ol.strip()), NEWCMD.match(nl_.strip())
                if m1 and m2:
                    rec.update(macro=m1.group(1), repo_value=m1.group(2), rev2_value=m2.group(2))
                else:
                    rec.update(kind="caption/note wording" if not rec["numbers_changed"] else "numbers",
                               repo_text=ol[:160], rev2_text=nl_[:160])
                changed.append(rec)
    # recompute the per-evaluation times from the run records
    sys.path.insert(0, HERE)
    import pandas as pd
    import mpce_results as R
    ALL, avail, fallbacks, new = R.load(HERE, False)
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(R.MAIN)].reset_index(drop=True)
    ms = G.Seconds / G.Calls * 1000
    med = ms.groupby(G.Algorithm).median()
    mean = ms.groupby(G.Algorithm).mean()
    tot = (G.Seconds.groupby(G.Algorithm).sum() / G.Calls.groupby(G.Algorithm).sum() * 1000)
    calls = G.Calls.groupby(G.Algorithm).agg(["min", "median", "max"])
    nruns = G.groupby("Algorithm").size()
    meta = [a for a in R.MAIN if a != "SLSQP" and a in med.index]
    per_n = ms.groupby([G.Algorithm, G.Turbines]).median().unstack()
    src = G.groupby("Algorithm").Source.agg(lambda s: ",".join(sorted(set(s))))
    sec_run = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack()
    cost = dict(
        data="6,030-evaluation runs, random initialization, data sets I and II (68 cases x 30 seeds per method); ms = Seconds / Calls * 1000 per run",
        runs=nruns.to_dict(), calls=calls.to_dict(orient="index"), source_files=src.to_dict(),
        median_ms={a: round(float(med[a]), 4) for a in med.index},
        mean_ms={a: round(float(mean[a]), 4) for a in mean.index},
        pooled_ms={a: round(float(tot[a]), 4) for a in tot.index},
        meta_median_range=[round(float(med[meta].min()), 4), round(float(med[meta].max()), 4)],
        meta_mean_range=[round(float(mean[meta].min()), 4), round(float(mean[meta].max()), 4)],
        slsqp_over_pso_median=round(float(med["SLSQP"] / med["PSOC"]), 4),
        slsqp_over_pso_mean=round(float(mean["SLSQP"] / mean["PSOC"]), 4),
        slsqp_over_pso_mean_from_rounded_table=round(round(float(mean["SLSQP"]), 3) / round(float(mean["PSOC"]), 3), 4),
        slsqp_over_pso_median_per_N={int(n): round(float(per_n.loc["SLSQP", n] / per_n.loc["PSOC", n]), 3) for n in per_n.columns},
        median_ms_per_N={a: {int(n): round(float(per_n.loc[a, n]), 3) for n in per_n.columns} for a in per_n.index},
        sec_per_run_range={a: [round(float(sec_run.loc[a].min()), 2), round(float(sec_run.loc[a].max()), 2)] for a in sec_run.index},
    )
    return changed, cost


INLINED = []


def inlined_tables():
    """Generated tables copied into the revision-2 section files (between "% >>>>> begin analysis/X.tex" and
    "% <<<<< end"): number sequence compared with the repository's generated file."""
    res = []
    for rel in MAIN_FILES:
        lines = open(os.path.join(R2, rel), encoding="utf-8").read().split("\n")
        i = 0
        while i < len(lines):
            m = re.match(r"%\s*>>>>> begin (analysis/\S+\.tex)", lines[i])
            if m:
                j = i + 1
                while j < len(lines) and not lines[j].startswith("% <<<<< end"):
                    j += 1
                blk = "\n".join(strip_comment(l) for l in lines[i + 1:j])
                rp = os.path.join(ROOT, m.group(1))
                ref = "\n".join(strip_comment(l) for l in open(rp, encoding="utf-8").read().split("\n"))
                a = [x["key"] if x["sign"] != "-" else "-" + x["key"] for x in extract_numbers(tex_to_plain(blk))]
                b = [x["key"] if x["sign"] != "-" else "-" + x["key"] for x in extract_numbers(tex_to_plain(ref))]
                diffs = []
                for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, b, a, autojunk=False).get_opcodes():
                    if tag != "equal":
                        diffs.append(dict(op=tag, repo=b[i1:i2], rev2=a[j1:j2]))
                res.append(dict(file=rel, lines=[i + 1, j + 1], source=m.group(1), n_numbers=len(a), n_repo=len(b),
                                identical_numbers=not diffs, diffs=diffs))
                i = j
            i += 1
    return res


def _table_by_label(text, label):
    i = text.find("\\label{%s}" % label)
    if i < 0:
        return None
    a = text.rfind("\\begin{table", 0, i)
    b = text.find("\\end{table", i)
    return "\n".join(strip_comment(l) for l in text[a:b].split("\n"))


def _nums(tex):
    return [("-" if x["sign"] == "-" else "") + x["key"] for x in extract_numbers(tex_to_plain(tex))]


REV2_TABLE_MAP = {"tab:S-rev-laplace": "tab:R2-laplace", "tab:S-rev-ga": "tab:R2-ga", "tab:S-rev-ga-hr": "tab:R2-ga-hr",
                  "tab:S-rev-ga-iea": "tab:R2-ga-iea", "tab:S-rev-spacing-rank": "tab:R2-spacing-rank",
                  "tab:S-rev-spacing": "tab:R2-spacing-contrast", "tab:S-rev-grad": "tab:R2-grad",
                  "tab:S-rev-grad-calls": "tab:R2-grad-calls", "tab:S-rev-site": "tab:R2-site"}
REV2_SECTIONS = [("sec:S-rev-lap-abl", "sec:S-rev-ga", ["/laplace"], ["tab:R2-laplace"]),
                 ("sec:S-rev-ga", "sec:S-rev-spacing", ["/ga68", "/ga_hr16", "/ga_iea37"], ["tab:R2-ga", "tab:R2-ga-hr", "tab:R2-ga-iea"]),
                 ("sec:S-rev-spacing", "sec:S-rev-grad", ["/spacing"], ["tab:R2-spacing-rank", "tab:R2-spacing-contrast"]),
                 ("sec:S-rev-grad", "sec:S-rev-site", ["/grad_iea37", "/ga_iea37"], ["tab:R2-grad", "tab:R2-grad-calls", "tab:R2-ga-iea"]),
                 ("sec:S-rev-site", "sec:S-further", ["/lillgrund"], ["tab:R2-site"])]


def rev2_tables_and_sections(rev2_out):
    """(a) The revision-2 tables typed into the supplement vs. the regenerated rev2_tables.tex (number sequences);
    (b) the prose of the revision-2 supplement sections against a section-scoped pool (the matching subtree of the
    regenerated rev2_summary.json and the section's regenerated tables): non-weak numbers not found are listed."""
    sup = open(os.path.join(R2, "SWEVO_supplement.tex"), encoding="utf-8").read()
    gen = open(os.path.join(rev2_out, "rev2_tables.tex"), encoding="utf-8").read()
    tabs = []
    for a, b in REV2_TABLE_MAP.items():
        x, y = _table_by_label(sup, a), _table_by_label(gen, b)
        na, nb = _nums(x or ""), _nums(y or "")
        diffs = [dict(op=t, regenerated=nb[i1:i2], supplement=na[j1:j2])
                 for t, i1, i2, j1, j2 in difflib.SequenceMatcher(None, nb, na, autojunk=False).get_opcodes() if t != "equal"]
        tabs.append(dict(supplement_table=a, regenerated_table=b, n_supplement=len(na), n_regenerated=len(nb),
                         identical_numbers=not diffs, diffs=diffs[:20]))
    summ = json.load(open(os.path.join(rev2_out, "rev2_summary.json")))
    flat = list(walk_json(summ))
    lines = sup.split("\n")
    lab_line = {m.group(1): i for i, l in enumerate(lines) for m in [re.search(r"\\label\{([^}]*)\}", l)] if m}
    prose = []
    for start, end, subtrees, tlabels in REV2_SECTIONS:
        pool = Pool()
        for p, x in flat:
            if any(p.startswith(st + "/") for st in subtrees):
                pool.add_float(x, p)
        for tl in tlabels:
            for k in _nums(_table_by_label(gen, tl) or ""):
                pool.add_key(k.lstrip("-"), tl)
        in_tab = False
        for i in range(lab_line[start], lab_line[end]):
            code = strip_comment(lines[i])
            if "\\begin{table" in code:
                in_tab = True
            if "\\end{table" in code:
                in_tab = False
                continue
            if in_tab or not code.strip():
                continue
            for n in extract_numbers(tex_to_plain(code)):
                if n["sig"] <= 2 and not n["exp"]:
                    continue
                k = n["key"]
                found = pool.get(k)
                if not found and n["pct"]:
                    found = pool.get(k)
                prose.append(dict(section=start, line=i + 1, text=("-" if n["sign"] == "-" else "") + k,
                                  found=bool(found), source=found[:2],
                                  context=re.sub(r"\s+", " ", tex_to_plain(code)[max(0, n["pos"] - 60):n["end"] + 40])))
    return tabs, prose


# ------------------------------------------------------------------------------------------------ part 2c (targeted)
def _rev2(exp):
    import pandas as pd
    fs = sorted(glob.glob(os.path.join(HERE, f"rev2_{exp}_s*of*.csv")))
    return pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in fs], ignore_index=True)


def _claim(out, where, claim, text, recomputed, ok, note=""):
    out.append(dict(where=where, claim=claim, text=text, recomputed=recomputed, status="OK" if ok else "DIFFERENT", note=note))


def targeted():
    """Recompute numbers of the revision-2 text that are not stored in any generated output (or only loosely)."""
    import numpy as np, pandas as pd
    out = []
    sup, beyond = "SWEVO_supplement.tex", "optA/sw/08_beyond.tex"
    # ---- elapsed-time table (tab:S-time) and the main-text ratios
    sp = _rev2("spacing")
    sp = sp[~((sp.Radius == 500) & (sp.Turbines == 10))]
    med_sp = sp.groupby("Algorithm").Seconds.median()
    lg = pd.concat([_rev2("lg16"), _rev2("lg16b")], ignore_index=True)
    med_lg = lg.groupby(["Budget", "Algorithm"]).Seconds.median()
    gr = _rev2("grad"); ga = _rev2("gaiea")
    gi = pd.concat([gr, ga], ignore_index=True)
    med_gi = gi.groupby(["Turbines", "Budget", "Algorithm"]).Seconds.median()
    rows = {}
    for a in ("PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "RSDVNS"):
        rows[a] = (f"{med_sp[a]:.2f} ({med_sp[a] / med_sp['PSOBV']:.2f})",
                   f"{med_lg[(6030, a)]:.2f} ({med_lg[(6030, a)] / med_lg[(6030, 'PSOBV')]:.2f})",
                   f"{med_lg[(30030, a)]:.2f} ({med_lg[(30030, a)] / med_lg[(30030, 'PSOBV')]:.2f})")
    typed = {"PSOBV": ("2.47 (1.00)", "5.68 (1.00)", "27.69 (1.00)"), "PSOC": ("2.07 (0.84)", "5.74 (1.01)", "28.94 (1.05)"),
             "SSABV": ("2.19 (0.89)", "5.56 (0.98)", "27.35 (0.99)"), "SSA": ("2.11 (0.85)", "5.70 (1.00)", "28.17 (1.02)"),
             "LXSSA": ("2.23 (0.90)", "5.69 (1.00)", "28.48 (1.03)"), "DE": ("2.90 (1.17)", "6.35 (1.12)", "31.18 (1.13)"),
             "BVNS": ("1.98 (0.80)", "5.58 (0.98)", "36.75 (1.33)"), "SLSQP": ("4.47 (1.81)", "12.79 (2.25)", "44.21 (1.60)"),
             "RSDVNS": ("2.14 (0.86)", "5.70 (1.00)", "31.67 (1.14)")}
    for a, t in typed.items():
        _claim(out, f"{sup}: tab:S-time", f"median s per run (ratio to PSO-VNS), {a}: spacing / Lillgrund 6,030 / 30,030",
               " | ".join(t), " | ".join(rows[a]), tuple(t) == rows[a],
               "spacing: 10 split cases (r=500, N=10 excluded), 5D and 6D pooled, 600 runs per method")
    for n, lab in ((16, "16"), (36, "36")):
        g6, g30 = med_gi[(n, 6030, "GA")], med_gi[(n, 30030, "GA")]
        for a, t in (("SLSQPX", {16: "4.3 / 31.8 (6.7)", 36: "15.4 / 84.8 (6.5)"}),
                     ("PSOSLSQPX", {16: "3.0 / 18.5 (3.9)", 36: "11.0 / 56.2 (4.3)"})):
            rec = f"{med_gi[(n, 6030, a)]:.1f} / {med_gi[(n, 30030, a)]:.1f} ({med_gi[(n, 30030, a)] / g30:.1f})"
            _claim(out, f"{sup}: tab:S-time", f"IEA37 {lab} turbines, {a}: median s 6,030 / 30,030 (ratio to GA at 30,030)",
                   t[n], rec, t[n] == rec)
        rec = f"{g6:.2f} / {g30:.2f} (1.00)"
        t = {16: "0.98 / 4.76 (1.00)", 36: "2.79 / 13.09 (1.00)"}[n]
        _claim(out, f"{sup}: tab:S-time", f"IEA37 {lab} turbines, GA median s", t, rec, t == rec)
    r_sp = med_sp["SLSQP"] / med_sp["PSOBV"]
    r6 = med_lg[(6030, "SLSQP")] / med_lg[(6030, "PSOBV")]; r30 = med_lg[(30030, "SLSQP")] / med_lg[(30030, "PSOBV")]
    _claim(out, f"{beyond}:271", "MS-SLSQP / PSO-VNS median time, spacing study", "1.8", f"{r_sp:.2f}", f"{r_sp:.1f}" == "1.8")
    _claim(out, f"{beyond}:271", "MS-SLSQP / PSO-VNS median time, Lillgrund 6,030 and 30,030", "2.3 and 1.6",
           f"{r6:.2f} and {r30:.2f}", (f"{r6:.1f}", f"{r30:.1f}") == ("2.3", "1.6"))
    rx = [med_gi[(n, 30030, "SLSQPX")] / med_gi[(n, 30030, "GA")] for n in (16, 36)]
    rp = [med_gi[(n, 30030, "PSOSLSQPX")] / med_gi[(n, 30030, "GA")] for n in (16, 36)]
    _claim(out, f"{beyond}:271", "IEA37 30,030: exact-gradient MS-SLSQP and PSO-SLSQP / GA median time", "6.5--6.7 and 3.9--4.3",
           f"{min(rx):.2f}--{max(rx):.2f} and {min(rp):.2f}--{max(rp):.2f}",
           (f"{min(rx):.1f}", f"{max(rx):.1f}", f"{min(rp):.1f}", f"{max(rp):.1f}") == ("6.5", "6.7", "3.9", "4.3"))
    # "Per run, the reported timings are about 3--6 times those of PSO-VNS" (exact-gradient methods, IEA37):
    # PSO-VNS IEA37 runs are in the repository (mpce_iea16p / mpce_iea36p)
    pv = pd.concat([pd.read_csv(os.path.join(HERE, f), usecols=lambda c: c not in ("Coordinates", "Curve"))
                    for f in ("mpce_iea16p_s0of1.csv", "mpce_iea36p_s0of1.csv")], ignore_index=True)
    pv = pv[(pv.Algorithm == "PSOBV") & (pv.Init == "random")]
    medp = pv.groupby(["Turbines", "Budget"]).Seconds.median()
    ratios = {f"{a} {n}T {b}": round(float(med_gi[(n, b, a)] / medp[(n, b)]), 2)
              for a in ("SLSQPX", "PSOSLSQPX") for n in (16, 36) for b in (6030, 30030)}
    mean_gi = gi.groupby(["Turbines", "Budget", "Algorithm"]).Seconds.mean()
    meanp = pv.groupby(["Turbines", "Budget"]).Seconds.mean()
    ratios_m = {f"{a} {n}T {b}": round(float(mean_gi[(n, b, a)] / meanp[(n, b)]), 2)
                for a in ("SLSQPX", "PSOSLSQPX") for n in (16, 36) for b in (6030, 30030)}
    lo, hi = min(ratios_m.values()), max(ratios_m.values())
    _claim(out, f"{beyond}: gradient paragraph (~l. 200)", "exact-gradient timings per run vs. PSO-VNS (ratio of mean "
           "times, the basis of sup:775; PSO-VNS IEA37 records mpce_iea16p/36p of the repository)", "about 3--6 times",
           f"means: {lo:.1f}--{hi:.1f} {ratios_m}; medians: {min(ratios.values()):.1f}--{max(ratios.values()):.1f}",
           2.5 <= lo and hi <= 6.5, "different batches / machines; the supplement (tab:S-time text) says the PSO-VNS IEA37 "
                                  "records belong to the missing archive, but they are in analysis/mpce_iea16p_s0of1.csv "
                                  "and mpce_iea36p_s0of1.csv (and the supplement itself uses them at l. 775)")
    # ---- gradient accounting
    g6 = gr[gr.Budget == 6030]
    for a in ("SLSQPX", "PSOSLSQPX"):
        for n in (16, 36):
            x = gr[(gr.Algorithm == a) & (gr.Turbines == n)]
            ok = ((x.FunCalls + x.CG * x.GradCalls + x.Unused) == x.Budget).all()
            _claim(out, f"{sup}: rev-grad", f"{a} {n}T: FunCalls + c_g GradCalls + Unused == Budget (all runs)", "yes",
                   str(bool(ok)), bool(ok))
    sx = g6[g6.Algorithm == "SLSQPX"]
    fm = sx.groupby("Turbines").FunCalls.mean(); gm = sx.groupby("Turbines").GradCalls.mean()
    _claim(out, f"{sup}:771", "exact-gradient MS-SLSQP at 6,030: mean objective / gradient evaluations (16, 36 turbines)",
           "1,667 and 1,647; 1,454 and 1,461", f"{fm[16]:.1f} and {fm[36]:.1f}; {gm[16]:.1f} and {gm[36]:.1f}",
           (round(fm[16]), round(fm[36]), round(gm[16]), round(gm[36])) == (1667, 1647, 1454, 1461))
    mx = (sx.FunCalls + 19 * sx.GradCalls).groupby(sx.Turbines).max()
    mx20 = (sx.FunCalls + 20 * sx.GradCalls).groupby(sx.Turbines).max()
    _claim(out, f"{sup}:771", "max over the 6,030 exact-gradient MS-SLSQP runs of FunCalls + 19 GradCalls (16, 36)",
           "29,486 and 29,550 (< 30,030); fits for any c_g <= 19", f"{mx[16]} and {mx[36]} (c_g = 20: {mx20[16]} and {mx20[36]})",
           (int(mx[16]), int(mx[36])) == (29486, 29550) and mx.max() <= 30030)
    sh = (sx.CG * sx.GradCalls / sx.Budget * 100).groupby(sx.Turbines).mean()
    sh30 = gr[(gr.Algorithm == "SLSQPX") & (gr.Budget == 30030)]
    sh30 = (sh30.CG * sh30.GradCalls / sh30.Budget * 100).groupby(sh30.Turbines).mean()
    _claim(out, f"{sup}:771", "share of the budget used by gradients, exact-gradient MS-SLSQP", "72--73%",
           f"6,030: {sh[16]:.1f} / {sh[36]:.1f}; 30,030: {sh30[16]:.1f} / {sh30[36]:.1f}",
           all(71.5 <= v < 73.5 for v in list(sh) + list(sh30)))
    # ---- SD ranges (sup:773)
    sds = gi[gi.Feasible.astype(str).str.lower().isin(["true", "1"])].groupby(["Algorithm", "Turbines", "Budget"]).Objective.std()
    gsd = [sds[(a, n, b)] for a in ("SLSQPX", "PSOSLSQPX") for n in (16, 36) for b in (6030, 30030)]
    pvf = pv[pv.Feasible.astype(str).str.lower().isin(["true", "1"])].groupby(["Turbines", "Budget"]).Objective.std()
    _claim(out, f"{sup}:773", "SD over seeds: exact-gradient methods vs. PSO-VNS (MWh)", "1,801--6,101 against 4,587--19,268",
           f"{min(gsd):.0f}--{max(gsd):.0f} against {pvf.min():.0f}--{pvf.max():.0f}",
           (round(min(gsd)), round(max(gsd)), round(pvf.min()), round(pvf.max())) == (1801, 6101, 4587, 19268))
    # ---- all-run paired outcome in the spacing study (tab:S-allrun-spacing)
    sp["F"] = sp.Feasible.astype(str).str.lower().isin(["true", "1"])
    key = ["Dataset", "Radius", "Turbines", "Seed", "Spacing"]
    typed_ar = {"PSOC": ("100.0", "178/0/122", "0.593", "80.0", "169/30/101", "0.613"),
                "SSABV": ("100.0", "236/0/64", "0.787", "80.7", "185/22/93", "0.653"),
                "SSA": ("91.3", "285/0/15", "0.950", "58.7", "251/28/21", "0.883"),
                "LXSSA": ("79.3", "293/0/7", "0.977", "52.0", "261/30/9", "0.920"),
                "DE": ("23.3", "294/0/6", "0.980", "10.0", "264/30/6", "0.930"),
                "BVNS": ("99.3", "251/0/49", "0.837", "82.0", "186/24/90", "0.660"),
                "SLSQP": ("100.0", "238/0/62", "0.793", "95.0", "171/6/123", "0.580"),
                "RSDVNS": ("100.0", "232/0/68", "0.773", "78.0", "212/20/68", "0.740"),
                ("SSABV", "RSDVNS"): ("--", "144/0/156", "0.480", "--", "143/44/113", "0.550")}
    for comp, t in typed_ar.items():
        first, second = (comp if isinstance(comp, tuple) else ("PSOBV", comp))
        rec = []
        for spc in ("5D", "6D"):
            A = sp[(sp.Algorithm == first) & (sp.Spacing == spc)].set_index(key)
            B = sp[(sp.Algorithm == second) & (sp.Spacing == spc)].set_index(key)
            j = A.join(B, lsuffix="_a", rsuffix="_b", how="inner")
            w = ((j.F_a & ~j.F_b) | (j.F_a & j.F_b & (j.Objective_a > j.Objective_b + 1e-9))).sum()
            l_ = ((~j.F_a & j.F_b) | (j.F_a & j.F_b & (j.Objective_b > j.Objective_a + 1e-9))).sum()
            tt = len(j) - w - l_
            feas = "--" if isinstance(comp, tuple) else f"{100 * B.F.mean():.1f}"
            rec += [feas, f"{w}/{tt}/{l_}", f"{(w + tt / 2) / len(j):.3f}"]
        _claim(out, f"{sup}: tab:S-allrun-spacing", f"{first} vs {second}: Feas. / W/T/L / score at 5D and 6D",
               " | ".join(t), " | ".join(rec), tuple(rec) == t, "bootstrap CIs not recomputed (seed not stored)")
    # ---- Lillgrund numbers from the records
    lgf = lg[lg.Feasible.astype(str).str.lower().isin(["true", "1"])]
    b30 = lgf[(lgf.Budget == 30030) & (lgf.Algorithm == "PSOBV")]
    import rev2_site_model as LGM
    inst = LGM.aep_gwh(LGM.site(16)[0])
    mxo = b30.Objective.max()
    _claim(out, f"{sup}: Lillgrund", "best PSO-VNS run at 30,030 exceeds installed block: AEP, % more",
           "116.16 GWh/yr, 0.05% more", f"{mxo:.2f}, {100 * (mxo / inst - 1):.2f}% (installed {inst:.2f}); runs above installed: "
           f"6,030 {int((lgf[lgf.Budget == 6030].Objective > inst).sum())}, 30,030 {int((lgf[lgf.Budget == 30030].Objective > inst).sum())}",
           f"{mxo:.2f}" == "116.16" and f"{100 * (mxo / inst - 1):.2f}" == "0.05")
    ms = b30.MinSpacing / 93.0
    _claim(out, f"{sup}: Lillgrund", "minimum spacing of the 30 PSO-VNS layouts at 30,030", "3.00D--3.08D",
           f"{ms.min():.2f}D--{ms.max():.2f}D (feasible: {len(b30)})", (f"{ms.min():.2f}", f"{ms.max():.2f}") == ("3.00", "3.08"))
    # ---- mean wall-clock times per run on IEA37 (sup:775)
    fdm = pd.concat([pd.read_csv(os.path.join(HERE, f), usecols=lambda c: c not in ("Coordinates", "Curve"))
                     for f in ("mpce_iea16_s0of1.csv", "mpce_iea36_s0of1.csv")], ignore_index=True)
    fdm = fdm[(fdm.Algorithm == "SLSQP") & (fdm.Init == "random")]
    order = [(16, 6030), (16, 30030), (36, 6030), (36, 30030)]
    mg = gr.groupby(["Algorithm", "Turbines", "Budget"]).Seconds.mean()
    mfd = fdm.groupby(["Turbines", "Budget"]).Seconds.mean()
    mpv = pv.groupby(["Turbines", "Budget"]).Seconds.mean()
    for lab, ser, t in (("exact-gradient MS-SLSQP", lambda k: mg[("SLSQPX",) + k], "4.3, 34.7, 16.0, 85.4"),
                        ("PSO-SLSQP", lambda k: mg[("PSOSLSQPX",) + k], "5.8, 23.7, 12.4, 56.6"),
                        ("MS-SLSQP forward differences", lambda k: mfd[k], "11.2, 56.8, 15.7, 78.8"),
                        ("PSO-VNS", lambda k: mpv[k], "1.3, 5.7, 3.2, 17.1")):
        rec = ", ".join(f"{ser(k):.1f}" for k in order)
        _claim(out, f"{sup}:775", f"mean s per run on IEA37 (16T 6,030 / 30,030; 36T 6,030 / 30,030), {lab}", t, rec, rec == t,
               "repository records mpce_iea16/36(p) for MS-SLSQP (FD) and PSO-VNS")
    # ---- IEA37 gaps to the best published layouts (strict 418,924.4 / 863,676.3; projected 421,451.1 / 882,382.8)
    pub = {("s", 16): 418924.4, ("s", 36): 863676.3, ("p", 16): 421451.1, ("p", 36): 882382.8}
    gg = lambda a, k: 100 * (1 - a / pub[k])
    _claim(out, f"{sup}:616 / {sup}:773", "GA best 36T at 30,030 (836,526.4): gap strict / projected", "3.14% / 5.20%",
           f"{gg(836526.4, ('s', 36)):.2f}% / {gg(836526.4, ('p', 36)):.2f}%",
           (f"{gg(836526.4, ('s', 36)):.2f}", f"{gg(836526.4, ('p', 36)):.2f}") == ("3.14", "5.20"))
    rec = (f"{gg(418757.4, ('s', 16)):.2f}", f"{gg(847780.3, ('s', 36)):.2f}", f"{gg(418757.4, ('p', 16)):.2f}", f"{gg(847780.3, ('p', 36)):.2f}")
    _claim(out, f"{sup}:773 / {beyond}", "best PSO-SLSQP 30,030 (418,757.4; 847,780.3): gap strict 16/36, projected 16/36",
           "0.04 / 1.84; 0.64 / 3.92", " / ".join(rec), rec == ("0.04", "1.84", "0.64", "3.92"))
    # ---- Lillgrund: jointly feasible PSO-VNS / PSO seeds at 30,030 (sup:849) and McNemar p values
    for b, tp in ((6030, "1.5e-5"), (30030, "9.8e-4")):
        A_ = lg[(lg.Budget == b) & (lg.Algorithm == "PSOBV")].set_index("Seed")
        B_ = lg[(lg.Budget == b) & (lg.Algorithm == "PSOC")].set_index("Seed")
        fa = A_.Feasible.astype(str).str.lower().isin(["true", "1"]); fb = B_.Feasible.astype(str).str.lower().isin(["true", "1"])
        n10, n01 = int((fa & ~fb).sum()), int((~fa & fb).sum())
        pm = min(1.0, 2 * sum(math.comb(n10 + n01, k) for k in range(min(n10, n01) + 1)) / 2 ** (n10 + n01))
        _claim(out, f"{sup}:849", f"Lillgrund {b}: exact McNemar p (PSO-VNS vs PSO feasibility)", tp,
               f"{pm:.1e} (discordant {n10}/{n01})", f"{pm:.1e}".replace("e-0", "e-") == tp)
        if b == 30030:
            j = fa & fb
            hi = int((A_.Objective[j] > B_.Objective[j]).sum())
            _claim(out, f"{sup}:849", "Lillgrund 30,030, jointly feasible seeds: n, PSO-VNS higher, means",
                   "19; 17; 114.60 vs 113.62", f"{int(j.sum())}; {hi}; {A_.Objective[j].mean():.2f} vs {B_.Objective[j].mean():.2f}",
                   (int(j.sum()), hi, f"{A_.Objective[j].mean():.2f}", f"{B_.Objective[j].mean():.2f}") == (19, 17, "114.60", "113.62"))
    # ---- Lillgrund parallelogram (sup:841): area, perimeter, Oler bound at 4D = 372 m
    res_g = []
    for fac, nm in ((1 / LGM.ENLARGE, "exact corner parallelogram"), (1.0, "boundary used (enlarged by 0.2%)")):
        poly = np.asarray(LGM.site(16)[1]) * fac
        xs, ys = poly[:, 0], poly[:, 1]
        area = 0.5 * abs(np.dot(xs, np.roll(ys, 1)) - np.dot(ys, np.roll(xs, 1)))
        per = np.sqrt((np.diff(np.r_[xs, xs[:1]]) ** 2 + np.diff(np.r_[ys, ys[:1]]) ** 2)).sum()
        d4 = 4 * 93.0
        oler = 2 * area / (np.sqrt(3) * d4 ** 2) + per / (2 * d4) + 1
        res_g.append((nm, f"{area / 1e6:.3f}", f"{per / 1e3:.2f}", f"{oler:.2f}"))
    _claim(out, f"{sup}:841", "Lillgrund parallelogram: area (km2), perimeter (km), Oler bound for 4D=372 m",
           "1.081; 4.24; 15.7", "; ".join(f"{a}: {b}, {c}, {d}" for a, b, c, d in res_g),
           any(r[1] == "1.081" for r in res_g), "conclusion (< 16 turbines) holds for both polygons")
    # ---- compact main-text equivalence table (tab:equiv-main) vs. the generated tab:X-equiv-levels and Bayes macros
    eqm = _table_by_label(open(os.path.join(R2, "optA/sw/06_results.tex"), encoding="utf-8").read(), "tab:equiv-main")
    lev = _table_by_label(open(os.path.join(HERE, "mpce_supp_inference.tex"), encoding="utf-8").read(), "tab:X-equiv-levels")
    MR_ = load_macros(sorted(glob.glob(os.path.join(HERE, "mpce_numbers*.tex"))))
    rope = {"PSO-VNS vs.\\ PSO": "NBayPSOVNSvsPSORope", "SSA-VNS vs.\\ RSD-VNS": "NBaySSAVNSvsRSDVNSRope",
            "SSA-VNS vs.\\ RS-VNS": "NBaySSAVNSvsRSVNSRope", "RSD-VNS vs.\\ RS-VNS": "NBayRSDVNSvsRSVNSRope",
            "LX-SSA-VNS vs.\\ RSD-VNS": "NBayLXSSAVNSvsRSDVNSRope", "PSO-VNS vs.\\ RSD-VNS": "NBayPSOVNSvsRSDVNSRope"}
    for row in eqm.split("\n"):
        pair = row.split("&")[0].strip()
        if pair not in rope:
            continue
        ref = [r for r in lev.split("\n") if r.split("&")[0].strip() == pair][0].split("&")
        # generated: n, dL, seed CI, m, Eq, case CI, m, p, Eq, CR2 CI, ...
        want = _nums("&".join([ref[2], ref[3], ref[6], ref[10]])) + _nums(expand("\\" + rope[pair], MR_))
        got = _nums(row)
        _claim(out, "optA/sw/06_results.tex: tab:equiv-main", f"{pair}: dL, seed / case / cluster 90% CI, P_rope",
               " ".join(got), " ".join(want), got == want, "tab:X-equiv-levels (analysis/mpce_supp_inference.tex) and \\NBay...Rope")
    # ---- revision record counts (design table and Section S-archive)
    cnt = {e: len(_rev2(e)) for e in ("ga", "gahr", "gaiea", "grad", "laplace", "spacing", "lg16", "lg16b")}
    _claim(out, f"{sup}: S-archive / tab:design", "revision record counts", "GA 2,040, GA HR 60, GA IEA37 120, grad 240, "
           "Laplace 2,160, spacing 6,480, Lillgrund 270 + 270 (11,640)", f"{cnt} total {sum(cnt.values())}",
           sum(cnt.values()) == 11640 and cnt["spacing"] == 6480 and cnt["laplace"] == 2160)
    # ---- equal-cluster table (tab:S-eqclus): case-weighted and equal-cluster means from the per-run data
    try:
        from scipy.stats import t as t_dist
        import mpce_inference_extra as MX
        G = MX.load(HERE)
        G = G[G.Feasible]
        cm = G.groupby(["Algorithm"] + MX.CASE).LossPct.mean()
        nf = G.groupby(["Algorithm"] + MX.CASE).size()
        typed_cl = {("PSOBV", "PSOC"): "0.064, 0.043, -0.008, -0.063, -0.031, -0.090",
                    ("SSABV", "RSDVNS"): "-0.153, 0.016, 0.025, -0.013, 0.064, 0.137",
                    ("SSABV", "RSVNS"): "-0.212, -0.055, -0.018, -0.119, -0.064, -0.008",
                    ("RSDVNS", "RSVNS"): "-0.059, -0.071, -0.042, -0.106, -0.128, -0.145",
                    ("LXBV", "RSDVNS"): "-0.067, 0.041, 0.059, 0.053, 0.152, 0.222",
                    ("PSOBV", "RSDVNS"): "-0.621, -0.233, -0.131, -0.454, -0.315, -0.235"}
        typed_eq = {("PSOBV", "PSOC"): ("-0.018", "-0.014", "[-0.063, 0.035]"),
                    ("SSABV", "RSDVNS"): ("+0.024", "+0.012", "[-0.067, 0.092]"),
                    ("SSABV", "RSVNS"): ("-0.068", "-0.079", "[-0.142, -0.017]"),
                    ("RSDVNS", "RSVNS"): ("-0.093", "-0.092", "[-0.125, -0.058]"),
                    ("LXBV", "RSDVNS"): ("+0.087", "+0.077", "[-0.005, 0.158]"),
                    ("PSOBV", "RSDVNS"): ("-0.306", "-0.332", "[-0.478, -0.185]")}
        for (a, b), t in typed_eq.items():
            d = (cm[a] - cm[b])
            q = (nf[a].reindex(d.index).fillna(0) >= 15) & (nf[b].reindex(d.index).fillna(0) >= 15)
            d = d[q]
            clm = d.groupby(level=[0, 1]).mean()
            m_eq = clm.mean(); se = clm.std(ddof=1) / np.sqrt(len(clm)); tc = t_dist.ppf(0.95, len(clm) - 1)
            rec = (f"{d.mean():+.3f}".replace("-0.", "-0.") if a != "PSOBV" or b != "PSOC" else f"{d.mean():.3f}",
                   f"{m_eq:+.3f}" if t[1].startswith("+") else f"{m_eq:.3f}", f"[{m_eq - tc * se:.3f}, {m_eq + tc * se:.3f}]")
            _claim(out, f"{sup}: tab:S-eqclus", f"{a} - {b}: case-weighted mean, equal-cluster mean, t5 90% CI",
                   " | ".join(t), " | ".join(rec), tuple(x.replace("+", "") for x in rec) == tuple(x.replace("+", "") for x in t),
                   f"n cases = {len(d)}")
            rc = ", ".join(f"{v:.3f}" for v in clm.values)
            _claim(out, f"{sup}: tab:S-eqclus", f"{a} - {b}: six cluster means (DS I 500/750/1000, DS II 500/750/1000)",
                   typed_cl[(a, b)], rc + "  [unrounded: " + ", ".join(f"{v:.5f}" for v in clm.values) + "]",
                   rc == typed_cl[(a, b)])
    except Exception as e:                                   # pragma: no cover
        _claim(out, f"{sup}: tab:S-eqclus", "equal-cluster table", "", f"not recomputed: {e!r}", False)
    return out


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rev2-out", default=None, help="directory with (or for) the regenerated rev2_summary.json / "
                    "rev2_tables.tex; default: a new temporary directory")
    ap.add_argument("--json", default=os.path.join(HERE, "rev3_numbers_audit.json"))
    a = ap.parse_args(argv)
    if a.rev2_out is None:
        import tempfile
        a.rev2_out = tempfile.mkdtemp(prefix="rev3_numbers_rev2out_")
    if not os.path.exists(os.path.join(a.rev2_out, "rev2_summary.json")):
        os.makedirs(a.rev2_out, exist_ok=True)
        subprocess.run([sys.executable, os.path.join(HERE, "rev2_analysis.py"), "--out-dir", a.rev2_out], check=True)
    M_repo = load_macros(sorted(glob.glob(os.path.join(HERE, "mpce_numbers*.tex"))))
    M_rev2 = load_macros(sorted(glob.glob(os.path.join(R2, "analysis", "mpce_numbers*.tex"))))
    stub = load_macros([os.path.join(R2, "optA", "phase4_macros_stub.tex"), os.path.join(R2, "optA", "phase6_macros_stub.tex")])
    print(f"macros: repository {len(M_repo)}, revision 2 {len(M_rev2)}, stubs {len(stub)} "
          f"(stub-only: {sorted(set(stub) - set(M_repo))})")
    macro_diff = {k: dict(repo=M_repo.get(k), rev2=M_rev2.get(k)) for k in sorted(set(M_repo) | set(M_rev2))
                  if M_repo.get(k) != M_rev2.get(k)}

    p1 = part1(M_repo, M_rev2)
    c1 = Counter(r["status"] for r in p1)
    print("part 1:", dict(c1))

    pools = build_pools(a.rev2_out)
    INLINED[:] = inlined_tables()
    p2_align, p2_nums = {}, []
    for rel in MAIN_FILES + SUPP_FILES:
        al = align(rel, M_repo, M_rev2)
        if al is not None:
            p2_align[rel] = al[0]
        p2_nums += scan_numbers(rel, M_rev2, pools, al)
    ca = Counter(r["status"] for v in p2_align.values() for r in v)
    print("part 2 alignment:", dict(ca))
    part_of = lambda r: "main" if r["file"] in MAIN_FILES else "supplement"
    cn = Counter((part_of(r), r["cls"]) for r in p2_nums)
    print("part 2 numbers:", dict(cn))

    p2t = INLINED
    print("part 2 inlined tables:", [(r["source"], r["identical_numbers"]) for r in p2t])
    r2tabs, r2prose = rev2_tables_and_sections(a.rev2_out)
    print("part 2 rev2 tables:", [(r["supplement_table"], r["identical_numbers"]) for r in r2tabs])
    print("part 2 rev2 prose:", Counter(r["found"] for r in r2prose))
    p2c = targeted()
    print("part 2c targeted:", Counter(r["status"] for r in p2c))
    p3, cost = part3()
    print("part 3 per-evaluation times:", json.dumps({k: cost[k] for k in ("median_ms", "mean_ms", "slsqp_over_pso_median",
                                                                           "slsqp_over_pso_mean")}, indent=0))
    out = dict(
        generated_by="analysis/rev3_numbers_audit.py",
        inputs=dict(repo_macros="analysis/mpce_numbers*.tex", rev2_macros="SWEVO_rev2/latex_source/analysis/mpce_numbers*.tex",
                    rev2_out=a.rev2_out, main_files=MAIN_FILES, supplement_files=SUPP_FILES),
        macro_files_differences=macro_diff,
        part1_check_comments=dict(counts=dict(c1), items=p1),
        part2_alignment=dict(counts=dict(ca), items={k: v for k, v in p2_align.items()}),
        part2_numbers=dict(counts={f"{k[0]}:{k[1]}": v for k, v in sorted(cn.items())}, items=p2_nums),
        part2_inlined_generated_tables=p2t,
        part2_rev2_tables_vs_regenerated=r2tabs,
        part2_rev2_supplement_prose=r2prose,
        part2c_targeted_recomputations=p2c,
        part3_hand_edits=dict(changed_lines=p3, cost_per_eval_recomputed=cost),
    )
    json.dump(out, open(a.json, "w"), indent=1, default=str)
    print("wrote", a.json)


if __name__ == "__main__":
    main()
