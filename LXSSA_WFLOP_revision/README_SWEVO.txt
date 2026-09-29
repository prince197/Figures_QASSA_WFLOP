SWEVO submission package (Swarm and Evolutionary Computation, Elsevier)
=====================================================================

Manuscript: "Baseline Configuration and Random-Sampling Controls in Metaheuristic Comparisons for Wind Farm
Layout Optimization" (P. Solanki, P. Dwivedi, V. Garg, V. Shukla).
This is a new submission. The MPCE version (MPCE_PSO_VNS.tex, MPCE_PSO_VNS_supplement.tex,
Cover_letter_PSO_VNS.tex) stays in the folder unchanged and still compiles (12 pages); both versions share the
number macros, generated tables and figures in analysis/ and figures_*/.
The SWEVO manuscript is the extended version (40 pages including references, 12-pt review format): it has its
own section files optA/sw/*.tex (Related Work; wind model and model figure; algorithm boxes; stability
proposition; parameter and experimental-design tables; definition of practical equivalence at three levels;
results with the coefficient sweep, equivalence levels and budget figure; feasibility lemma; recommendations
table; threats to validity). Detailed tables stay in the supplement; the log of what moved between the two
documents is optA/sw/moved.md (items moved back to the supplement to fit 40 pages are marked
"% Page budget" in optA/sw/*.tex). The reference list is ordered by first citation (analysis/swevo_bib.py).
Journal requirements (checklist, with what is verified and what is assumed): optA/swevo_requirements.md.


1. CONTENTS OF THE SUBMISSION ZIP
---------------------------------
Compiled documents (root):
  SWEVO_manuscript.tex / .pdf     main manuscript (elsarticle, preprint 12pt, line numbers, numbered refs)
  SWEVO_supplement.tex / .pdf     supplementary material (Tables S*, Figs. S*, theory)
  SWEVO_cover_letter.tex / .pdf   cover letter to the Editors-in-Chief
  SWEVO_declarations.tex / .pdf   convenience draft: competing interest, CRediT, funding, data, generative AI
  SWEVO_highlights.txt            3-5 highlights, one per line, each <= 85 characters incl. spaces (plain text)
  SWEVO_manuscript_full.tex       the SAME main manuscript as ONE LaTeX file (every \input, number macro and
                                  generated table inlined; analysis/flatten_swevo.py). It needs only figures_mpce/,
                                  figures_final/ and SWEVO_supplement.aux (or xref/, see 2b) and gives the identical
                                  40-page PDF. Regenerate it after any edit of the modular files.
Sources read by the .tex files:
  optA/swevo_front.tex            title, authors, affiliations, abstract, keywords, nomenclature (SWEVO only)
  optA/sw/*.tex                   sections of the SWEVO version (02_intro ... 10_limits_concl; capture.tex
                                  places generated supplement tables/figures in the main text)
  optA/swevo_back.tex             declarations (CRediT, competing interest, data, generative AI, funding/
                                  acknowledgments) and the reference list (SWEVO only)
  optA/supp_theory.tex            theory section of the supplement
  optA/phase4_macros_stub.tex, optA/phase6_macros_stub.tex   additional number macros
  analysis/mpce_numbers.tex, _extra, _diag, _dir, _csweep   number macros (generated; do not edit by hand)
  analysis/mpce_tab_*.tex, analysis/mpce_supp_*.tex, analysis/mpce_supplementary.tex   generated tables
  figures_mpce/*.pdf, figures_final/*.pdf                    figures
  xref/ (Overleaf / upload copy only; see 2b)               frozen .aux files for the cross-references
Not compiled, for reference: optA/SWEVO.md (brief), optA/swevo_requirements.md, this README.
Do not include in the ZIP: the MPCE files, Response_to_previous_reviews_PSO_VNS.md, *.log, *.aux (except xref/).


2. HOW TO COMPILE
-----------------
Compiler: pdfLaTeX (elsarticle.cls is part of TeX Live and of Overleaf). No BibTeX run is needed: the
reference list is a thebibliography inside optA/swevo_back.tex.

The manuscript cites Tables S*/Figs. S* of the supplement and the supplement cites sections, tables and
equations of the manuscript, through the package xr (\externaldocument). Each document reads the .aux file of
the other, so they must be compiled alternately.

2a. Local build (recommended): from this folder run
      sh analysis/build_swevo.sh
    It compiles SWEVO_supplement.tex and SWEVO_manuscript.tex alternately (three rounds), then the cover
    letter and the declarations, and prints pages, undefined references and overfull boxes for each PDF.
    By hand:  pdflatex SWEVO_supplement; pdflatex SWEVO_manuscript; repeat both twice more;
              pdflatex SWEVO_cover_letter (twice); pdflatex SWEVO_declarations (twice).
    The build does not rerun the Python analysis. To regenerate numbers, tables and figures, run
    sh analysis/build_mpce_paper.sh or the scripts listed in analysis/README_reproduce.md first.

2b. Overleaf (same xref/ trick as README_OVERLEAF.txt of the MPCE project): Overleaf compiles one main
    document at a time and cannot read the .aux file of the other document directly. Therefore
      - create a folder xref/ and put in it the current SWEVO_supplement.aux and SWEVO_manuscript.aux
        (from a local build);
      - in the Overleaf copy only, change the two xr lines to
          SWEVO_manuscript.tex:   \externaldocument[S-]{xref/SWEVO_supplement}
          SWEVO_supplement.tex:   \externaldocument[M-]{xref/SWEVO_manuscript}
      - Menu > Compiler: pdfLaTeX; Menu > Main document: SWEVO_manuscript.tex (or the supplement, cover
        letter or declarations) > Recompile.
    If you change the numbering of sections, tables, figures or equations in one document, compile it,
    download its .aux (Logs and output files > Other logs and files > .aux), replace the file of the same
    name in xref/ and recompile the other document.
    The same xref/ copy is the one to put in a LaTeX source ZIP for Editorial Manager, whose PDF builder also
    compiles one file at a time.

After compiling, check that no "??" appears, and that the only red [TBD: ...] text is the author
information of Section 4.


3. WHAT TO UPLOAD IN ELSEVIER EDITORIAL MANAGER
-----------------------------------------------
Article type: research article ("Research paper" / "Full Length Article"; exact label: see EM).
  Item type in EM               File
  Manuscript                    SWEVO_manuscript.pdf, plus the LaTeX source ZIP (Section 1, with xref/ as in
                                2b; EM may list it as "LaTeX source files"). If EM accepts a single PDF at
                                first submission ("Your Paper Your Way"), the source is required at revision.
  Supplementary material        SWEVO_supplement.pdf (caption, e.g., "Supplementary material: model details,
                                per-case results, further tables and figures, and theory")
  Highlights                    SWEVO_highlights.txt (editable file; "Highlights" in the file name)
  Cover letter                  SWEVO_cover_letter.pdf (or paste its text into the cover-letter field)
  Declaration of interest       REQUIRED: a .docx generated by the authors with Elsevier's Declaration of
                                Interests tool ("I have nothing to declare" if that applies). It cannot be
                                produced from this folder.
  Declarations (optional)       SWEVO_declarations.pdf is a convenience draft collecting the competing-interest,
                                CRediT, funding, data and generative-AI statements (all also in the manuscript
                                back matter); use it to fill in the EM forms, or upload it if EM offers a slot.
No separate title page and no anonymized manuscript: review at SWEVO is single-anonymized, so authors and
affiliations stay in SWEVO_manuscript.pdf.
In the EM forms also enter: title, abstract and keywords (copy from optA/swevo_front.tex), the funding
information, the data availability statement (Zenodo DOI), ORCID of all authors (at least the corresponding
author), and, if wanted, suggested reviewers. A graphical abstract is optional (not prepared).


4. REMAINING AUTHOR ITEMS (red [TBD: ...] text)
-----------------------------------------------
  - Corresponding author: confirm (optA/swevo_front.tex, cover letter, declarations).
  - Funding statement (optA/swevo_back.tex, SWEVO_declarations.tex); if there is none, Elsevier's wording:
    "This research did not receive any specific grant from funding agencies in the public, commercial, or
    not-for-profit sectors."
  - Acknowledgment: write or delete (optA/swevo_back.tex).
  - CRediT roles: confirm or edit for all four authors (optA/swevo_back.tex and SWEVO_declarations.tex; keep
    both identical).
  - Declaration of competing interest: confirm "no competing interests" (manuscript, cover letter,
    declarations).
  - Generative-AI declaration: confirm or edit the wording (manuscript and declarations). Heading used:
    "Declaration of generative AI and AI-assisted technologies in the manuscript preparation process"
    (Elsevier's current wording).
  - Declaration of competing interest: generate the .docx with Elsevier's Declaration of Interests tool.
  - Zenodo DOI of code and data (optA/swevo_back.tex, cover letter, declarations).
  - Volume, issue and pages of Solanki & Deep 2023 (reference list, optA/swevo_back.tex).
  - Cover letter: date of submission; suggested reviewers (optional; delete the sentence if none). The
    salutation is "Dear Editors-in-Chief," without names, because the current Editors-in-Chief could not be
    verified (see optA/swevo_requirements.md); add names only after checking the editorial-board page.
  - Optional: a paragraph disclosing the earlier MPCE submission of a related manuscript, if the authors wish
    (the MPCE letter, Cover_letter_PSO_VNS.tex, has such a paragraph; it is not in the SWEVO letter).
  - Numbers whose runs are still pending, if any, print as [TBD: pending ...]; see the second line of
    analysis/mpce_numbers.tex and the "pending" count printed by analysis/build_swevo.sh.


5. TO CONFIRM ON THE GUIDE FOR AUTHORS BEFORE SUBMITTING
--------------------------------------------------------
The requirements were collected from search snippets only (the Elsevier pages were not reachable); see
optA/swevo_requirements.md and check at
https://www.sciencedirect.com/journal/swarm-and-evolutionary-computation/publish/guide-for-authors :
  - Abstract length limit (assumed about 250 words; check the abstract in optA/swevo_front.tex against it).
  - Research-data option (reportedly Option C: deposit, cite and link the data, or state why not) and the
    wording of the data availability statement.
  - Whether a CRediT authorship statement is required at submission or only at revision.
  - Whether suggested reviewers are required, and how many (EM usually asks for 3-5).
  - Article processing charge for open access (reportedly about USD 3,160 excl. taxes; hybrid journal) and
    whether an institutional agreement covers it; subscription publication has no charge.
  - Also: the exact article-type label in EM, the current Editors-in-Chief, and the current heading of the
    generative-AI declaration.
