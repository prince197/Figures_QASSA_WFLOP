# Conversion to Swarm and Evolutionary Computation (Elsevier) — brief for all agents

Goal: a submission-ready package for Swarm and Evolutionary Computation (SWEVO), keeping the MPCE version
(MPCE_PSO_VNS.tex, MPCE_PSO_VNS_supplement.tex) compiling and unchanged in content.

Files (new; each agent owns only the files named in its task):
- SWEVO_manuscript.tex — elsarticle main file: preamble (packages as in MPCE_PSO_VNS.tex minus IEEE items; number
  macros \input in the same order: analysis/mpce_numbers.tex, _extra, _diag, _dir, _csweep, optA/phase4_macros_stub.tex,
  optA/phase6_macros_stub.tex; \TBD and \pipefig definitions; xr: \externaldocument[S-]{SWEVO_supplement}),
  then \input{optA/swevo_front} + \input{optA/02_intro} … \input{optA/10_limits_concl} + \input{optA/swevo_back}.
- optA/swevo_front.tex — elsarticle frontmatter: title, authors/affiliations/emails (from optA/01_front.tex
  \thanks), corresponding author (\TBD), abstract (from 01_front, ≤ 250 words, no citations, no CHECK comments
  inside the abstract body are fine as comments), keywords (\begin{keyword} … \sep …, 1–7), highlights not in the
  manuscript (separate file), nomenclature as a compact table or description list (not IEEEdescription).
- optA/swevo_back.tex — CRediT authorship statement (\TBD roles), Declaration of competing interest, Data availability,
  Declaration of generative AI and AI-assisted technologies in the manuscript preparation process (Elsevier policy; wording given by
  the lead — see below), Acknowledgments/Funding (\TBD), references: the same \bibitem list as optA/11_back.tex
  (keep keys and order; numbered style [1]; no biographies).
- SWEVO_supplement.tex — copy of MPCE_PSO_VNS_supplement.tex retargeted (xr prefix M- → SWEVO_manuscript; no
  journal-specific wording).
- SWEVO_cover_letter.tex, SWEVO_highlights.txt, SWEVO_declarations.tex (competing interests + CRediT + AI
  statement as a separate form-like document), SWEVO_title_page.tex (only if double-anonymized review applies),
  README_SWEVO.txt.
- analysis/build_swevo.sh — builds supplement and manuscript alternately (xr) + cover letter + declarations.

Rules
- Section files optA/02_intro.tex … optA/10_limits_concl.tex are shared by both versions: edits must work in both
  (IEEEtran and elsarticle). Only the cleanup agent edits them.
- Numbers only via the existing macros. Keep all CHECK comments.
- Do not run pdflatex in the shared folder; compile privately in the scratchpad:
  /tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/sw_<you>/
- No commit/push. Report files, compile result (pages, errors, undefined refs, overfull boxes).

Generative-AI declaration (lead wording; authors must confirm; keep the \TBD marker):
"During the preparation of this work the authors used Claude (Anthropic) to assist with code for the analyses,
drafting and editing of the text, and consistency checks. After using this tool, the authors reviewed and edited
the content as needed and take full responsibility for the content of the published article.
\TBD{authors: confirm or edit this statement}"
