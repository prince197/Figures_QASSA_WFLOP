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
    prim, sec, macro = Pool(), Pool(), Pool()
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
    jsons = [(os.path.join(rev2_out, "rev2_summary.json"), "rev2_summary.json(regenerated)")]
    jsons += [(f, os.path.relpath(f, ROOT)) for f in sorted(glob.glob(os.path.join(VAL, "*.json")))
              if not f.endswith("raw_data_checksums.json")]
    jsons += [(f, os.path.relpath(f, ROOT)) for f in sorted(glob.glob(os.path.join(HERE, "mpce_summary*.json")))]
    for fn, lab in jsons:
        for p, x in walk_json(json.load(open(fn))):
            prim.add_float(x, f"{lab}{p}")
    for fn in sorted(glob.glob(os.path.join(VAL, "*.csv"))):
        for i, row in enumerate(csv.DictReader(open(fn)), 2):
            for k, v in row.items():
                prim.add_float(v, f"{os.path.relpath(fn, ROOT)}:{i}:{k}")
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
    return M, macro, prim, sec


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
    return results, nt, new_owner, nl


# ------------------------------------------------------------------------------------------------ part 2b (all numbers)
YEARS = set(range(1950, 2031))
STRUCT_INT = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 24, 25, 30, 36, 40, 45, 50, 60, 64, 68,
              80, 90, 100, 120, 150, 180, 200, 250, 300, 308, 360, 385, 400, 462, 500, 600, 750, 1000, 1300, 2000,
              3015, 6030, 30030, 120030, 3000, 6000, 30000, 120000, 1000000}
STRUCT_DEC = {"0.5", "0.05", "0.95", "0.9", "0.75", "0.25", "1.5", "2.5", "0.7", "0.7298", "1.49618", "1.496", "2.0",
              "0.01", "0.1", "0.2", "0.04", "0.075", "0.3", "0.6", "0.8", "0.4", "1.0", "0.0", "0.001"}


def scan_numbers(rel, M_rev2, pools, align_info):
    M_unused, macro, prim, sec = pools
    fn = os.path.join(R2, rel)
    text = open(fn, encoding="utf-8").read()
    out = []
    # per line plain text (macros expanded with the printed = rev2 values)
    nt, new_owner, nl = (align_info[1], align_info[2], align_info[3]) if align_info else (None, None, None)
    anchored = defaultdict(list)    # line -> [(token, macro)]
    if nt is not None:
        for t, o, l in zip(nt, new_owner, nl):
            if o and re.match(r"\d", t):
                anchored[l].append((t, o))
    for no, line in enumerate(text.split("\n"), 1):
        code = strip_comment(line)
        if not code.strip():
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
            if anc.get(k):
                anc[k] -= 1
                rec.update(cls="anchored_equal", source=["\\" + anc_m[k].pop(0)])
            elif used_macros and any(k in [x["key"] for x in extract_numbers(tex_to_plain(expand("\\" + u, M_rev2)))]
                                     for u in used_macros):
                rec.update(cls="macro_in_rev2_text", source=["\\" + u for u in used_macros
                                                              if k in [x["key"] for x in extract_numbers(tex_to_plain(expand("\\" + u, M_rev2)))]][:3])
            else:
                src = macro.get(k) + prim.get(k)
                if src:
                    rec.update(cls="traced", source=src[:4])
                elif sec.get(k):
                    rec.update(cls="traced_secondary", source=sec.get(k)[:3])
                else:
                    v = abs(n["value"])
                    struct = (not n["exp"] and n["dec"] == 0 and (v in STRUCT_INT or int(v) in YEARS or v <= 12)) \
                        or k in STRUCT_DEC
                    rec.update(cls="structural" if struct else "untraced", source=[])
                    if not struct:
                        rec["near"] = near_values(n, macro, prim)
            rec["weak"] = n["sig"] <= 2 and not n["exp"]
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


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rev2-out", default="/tmp/rev3_numbers_rev2out")
    ap.add_argument("--json", default=os.path.join(HERE, "rev3_numbers_audit.json"))
    a = ap.parse_args(argv)
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
        part3_hand_edits=dict(changed_lines=p3, cost_per_eval_recomputed=cost),
    )
    json.dump(out, open(a.json, "w"), indent=1, default=str)
    print("wrote", a.json)


if __name__ == "__main__":
    main()
