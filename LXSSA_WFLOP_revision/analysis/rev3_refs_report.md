# Revision 3: reference audit (A5, A6) and SWEVO Guide-for-Authors check (A7)

Date: 2026-10-04. Topic: `refs`. Files: this report and `analysis/rev3_refs_corrections.tex` (corrected `\bibitem`s only).
I did not edit any source.

## 0. How the checks were done and what they rest on

- The egress proxy blocked every primary bibliographic host I tried, for both `curl` and WebFetch (HTTP 403 /
  EGRESS_BLOCKED): doi.org, api.crossref.org, api.openalex.org, api.semanticscholar.org, dblp.org, link.springer.com,
  sciencedirect.com, elsevier.com, wes.copernicus.org, arxiv.org, jmlr.org, ieeexplore, wiley, sage, mit, mdpi,
  researchgate, ideas.repec.org and DTU Orbit. **Crossref was not reachable, so no entry was checked against Crossref.**
- Reachable sources:
  - **S**: WebSearch result records, i.e. search-engine summaries of publisher, ADS, OSTI, DTU Orbit, ResearchGate
    and similar pages. This is secondary evidence. Each S check named the fields in the query and compared the
    returned authors, volume, issue, pages, year and DOI with the entry.
  - **PyPI** JSON API (direct).
  - **git** protocol to github.com (`git ls-remote` and `git fetch`), used for the IEA37 repository, the JAX README
    and the Python, NumPy and SciPy release tags.
  - **ENV**: the local environment (`python3 --version`, `numpy`, `scipy` and `py_wake` versions).
- **P** means the field is covered by an earlier `% VERIFY-REF (2026-09-28/29/30)` comment that names primary
  records. I re-checked those entries for internal consistency but could not reopen the primary pages.
- **K** means the record agrees with the standard bibliographic record of a classic paper. I did not re-fetch it in
  this session; for most of these, an S check was also done.

Status codes:
- **OK**: all fields agree.
- **corrected**: a change is supplied in `rev3_refs_corrections.tex`.
- **unverifiable**: a field could not be confirmed. The field is named in the row.

## 1. A5: Solanki & Deep (2023)

- **Title.** The publisher's form is "Laplacian Salp Swarm Algorithm for continuous optimization". The bibliography
  uses sentence case, as for the other entries.
- **Authors.** Prince Solanki and Kusum Deep, Department of Mathematics, IIT Roorkee.
- **Journal and DOI.** International Journal of System Assurance Engineering and Management, doi
  10.1007/s13198-023-01935-y.
- **History.** Received 18 Aug 2021, revised 2 Feb 2023, accepted 1 May 2023, published online 18 May 2023.
- **Source.** Search-index records of the Springer article page and of the journal's "Online first articles"
  listing (page 8).
- **Volume, issue and pages: not found.** The indexed publisher citation reads "Int J Syst Assur Eng Manag (2023).
  https://doi.org/10.1007/s13198-023-01935-y" with no volume. The article still appears in the journal's Online
  First listing. Crossref and the Springer page were blocked, so I could not confirm this directly.
- **Correction.** I replaced the IEEE term "early access, May 2023" with "published online May 18, 2023 (Online
  First)". The supplement entry was made the same.
- **Action for the first author.** Check the Springer page at proof stage and add the volume, issue and pages if
  they have been assigned by then.

## 2. A6: audit of the main bibliography (`optA/swevo_back.tex`, 88 entries)

- All 88 keys are cited in the manuscript, and every citation has an entry (script check).

| # | Key | Status | What differs / note | Source |
|---|-----|--------|---------------------|--------|
| 1 | Srikakulapu2018 | OK | vol. 6, pp. 1181–1192, 2018 confirmed; no. 6 / Nov. not independently re-confirmed | S, K |
| 2 | Hou2019 | OK | 6 authors, 7(5) 975–986, DOI | S (DTU Orbit, Springer) |
| 3 | Kusiak2010 | OK | 35(3) 685–694, DOI | S (OSTI, RG) |
| 4 | Baker2019 | OK | 5 authors, AIAA 2019-0540, DOI 10.2514/6.2019-0540 | S |
| 5 | Azlan2021 | OK | — | P |
| 6 | Pookpunt2013 | OK | — | P |
| 7 | Asaah2021 | OK | 9(2) 367–375, 2021 | S |
| 8 | Mirjalili2017 | OK | — | K |
| 9 | Rezk2019 | OK | 12(22) art. 4335 | S (MDPI page) |
| 10 | Hooker1995 | OK | — | P |
| 11 | Derrac2011 | OK | — | P |
| 12 | BartzBeielstein2020 | OK | arXiv:2007.03488 (v2, 16 Dec 2020); no journal version found | S |
| 13 | LaTorre2021 | OK | — | P |
| 14 | Sorensen2015 | OK | 22(1) 3–18, DOI | S (Wiley) |
| 15 | Aranha2022 | OK | — | P |
| 16 | Wilson2018 | OK | vol. 126 (Oct. 2018) 681–691; first 10 authors listed, "et al." correct | S, P |
| 17 | Thomas2023 | **corrected** | "no. 5" removed. The publisher's citation "Wind Energ. Sci., 8, 865–891, 2023" has no issue. Secondary indexes give 8(5), but the article falls between 8/691 (5 May, issue 5) and 8/947 (7 Jun, issue 6), so the issue is unconfirmed. The old comment said "issue omitted" while the entry kept no. 5. Main and supplement. | S (Copernicus, OSTI, ADS, ProQuest) |
| 18 | Kennedy1995 | OK | ICNN'95 vol. 4 pp. 1942–1948, DOI | S |
| 19 | Poli2009 | OK | 13(4) 712–721 | S |
| 20 | Bonyadi2017 | OK | — | P |
| 21 | Cleghorn2018 | OK | vol. 12 pp. 1–22, DOI; no. 1 from K | S |
| 22 | Benavoli2016 | OK | JMLR 17(5):1–10 (no DOI) | S (jmlr.org) |
| 23 | Lakens2017 | OK | — | P |
| 24 | Mladenovic1997 | OK | — | K |
| 25 | Hansen2001 | OK | — | K |
| 26 | Solanki2023 | **corrected** (vol/no/pp **unverifiable**) | See Section 1 | S |
| 27 | Neri2012 | OK | — | P |
| 28 | Deb2000 | OK | 186(2–4) 311–338 | S |
| 29 | Wu2017 | OK | CEC 2017 constrained tech. report, 2017 | S |
| 30 | Benavoli2017 | OK | — | P |
| 31 | Mosetti1994 | OK | 51(1) 105–116 | S |
| 32 | Grady2005 | OK | 30(2) 259–270 | S |
| 33 | Jensen1983 | OK | Risø-M-2411 | K |
| 34 | Katic1986 | OK (note) | DTU Orbit dates the proceedings volume 1987 (EWEC'86, Rome, vol. 1, pp. 407–410). Conference date Oct. 1986 kept, per the IEEE convention. Optional: add "(publ. 1987)". | S (DTU Orbit) |
| 35 | Bansal2017 | OK | vol. 107, Jul. 2017, 386–402 | S |
| 36 | Eroglu2012 | OK | 44, 53–62 | S |
| 37 | Eroglu2013 | OK | 58, 95–107 | S |
| 38 | Feng2015 | OK | — | P |
| 39 | Elkinton2008 | OK | 32(1) 67–84 | S |
| 40 | Ju2019 | OK | 248, 429–445, DOI | S |
| 41 | Bai2022 | OK | 5 authors, Jan. 2022, art. 115047 (ADS 2022ECM...25215047B) | S |
| 42 | Nagpal2021 | OK | Shriya V. Nagpal, M. Vivienne Liu, C. Lindsay Anderson; 168, May 2021, 581–592 | S |
| 43 | Cazzaro2022 | OK | 138, art. 105588 | S |
| 44 | IEA37repo | OK (comment corrected) | Commit 267f6e5 is still HEAD of IEAWindSystems/iea37-wflo-casestudies (13 Nov 2020; merge of PR #1 from byuflowlab). Its tree is identical to byuflowlab af88908. The byuflowlab repository still exists separately, so the old comment's "redirects" is wrong; only the comment is refreshed. | git |
| 45 | Guirguis2016 | OK | — | P |
| 46 | Quaeghebeur2021 | OK | 6(3) 815–839 | S |
| 47 | Stanley2019 | OK | — | P |
| 48 | Bastankhah2014 | OK | 70, 116–123 | S |
| 49 | Tao2020 | OK | 5 authors, 159, Oct. 2020, 553–569 | S |
| 50 | Rodrigues2024 | OK | 9(2) 321–341, 5 authors | S |
| 51 | Quick2023 | OK | 8(8) 1235–1250, 5 authors | S |
| 52 | CriadoRisco2024 | OK | 6 authors, 9, 585–600; "no. 3" fits a Mar. 2024 publication but is not independently confirmed | S |
| 53 | Yang2023 | OK | 6 authors, 218, Dec. 2023, art. 119240 | S |
| 54 | Wang2024 | OK | 6 authors, 293, art. 116644 | S |
| 55 | Li2025 | OK | 6 authors, 392, Aug. 2025, art. 125908 | S |
| 56 | ParkPark2019 | OK | 187, art. 115883 | S |
| 57 | LiZhang2023 | OK | 339 (1 Jun 2023), art. 120928 | S |
| 58 | LiRobert2024 | OK | 359 (Apr. 2024), art. 122758 | S |
| 59 | Bempedelis2024 | OK | 9(4) 869–882 | S |
| 60 | Hansen2021 | OK | — | P |
| 61 | Demsar2006 | OK | — | K |
| 62 | Holm1979 | OK | — | K |
| 63 | Carrasco2020 | OK | — | P |
| 64 | Campelo2019 | OK | — | P |
| 65 | CamachoVillalon2023 | OK | — | P |
| 66 | Castelli2022 | OK | 189, art. 116029 | S |
| 67 | Shi1998 | OK | pp. 69–73, DOI | S |
| 68 | Clerc2002 | OK | 6(1) 58–73, DOI | S |
| 69 | Trelea2003 | OK | — | P |
| 70 | VandenBergh2006 | OK | — | P |
| 71 | Krasnogor2005 | OK | — | P |
| 72 | Mladenovic2008 | OK | 191(3) 753–770 | S |
| 73 | Tasgetiren2007 | OK | 177, 1930–1947 | S |
| 74 | Marinakis2017 | OK | 261(3) 819–834 | S |
| 75 | Gumaida2019 | OK | 49, 3539–3557 (online 18 Apr 2019); no. 10 from K | S |
| 76 | Gebraad2016 | OK | 19(1) 95–114 | S |
| 77 | Bastankhah2016 | OK | 806, 506–541 | S |
| 78 | Zaharie2002 | OK | MENDEL 2002, Brno, pp. 62–67 | S |
| 79 | Deep2007 | OK | 188(1) 895–911 | S |
| 80 | Storn1997 | OK | 11, 341–359; no. 4 from K | S |
| 81 | Kraft1988 | OK | — | K |
| 82 | Deb1995 | OK | 9(2) 115–148 | S |
| 83 | DebGoyal1996 | OK | 26(4) 30–45 | S |
| 84 | PyWake | OK | py-wake 2.6.20 uploaded 2026-03-30, MIT License, author "DTU Wind Energy", source URL as cited; 2.6.20 installed here | PyPI, ENV |
| 85 | Gocmen2016 | OK | — | P |
| 86 | Carrillo2013 | OK | 21, 572–581 | S |
| 87 | HansenOstermeier2001 | OK | 9(2) 159–195 | S |
| 88 | Tanabe2014 | OK | CEC 2014 pp. 1658–1665, DOI | S |

Summary: 85 OK (one with an optional note), 2 corrected (Thomas2023, Solanki2023) and 1 comment-only correction
(IEA37repo). One field set is unverifiable: the volume, issue and pages of Solanki2023, which have not been
assigned.

## 3. Supplement bibliography (`SWEVO_supplement.tex`, 23 entries)

| Key | Status | Note |
|-----|--------|------|
| Kusiak2010, Eroglu2012, Eroglu2013, Bansal2017, Mirjalili2017, Poli2009, Kraft1988, Bastankhah2014, Baker2019, IEA37repo, Benavoli2017, Lakens2017, Deb1995, DebGoyal1996, Guirguis2016, Gocmen2016 | OK | Same fields as the main entries. Some are abbreviated (no month, "et al."), which is acceptable. |
| Solanki2023 | **corrected** | Same text as the main entry |
| Thomas2023 | **corrected** | "no. 5" removed, as in the main entry |
| PyWake | **corrected (style)** | "Accessed" moved before "[Online]" to match the main entry |
| JAX2018 | OK | The README of tag jax-v0.10.2 (commit 990e6a0, 17 Jun 2026) gives the same 12 authors, title and year 2018. PyPI: jax 0.10.2 uploaded 2026-06-17. The archive has no file importing JAX, so the JAX cross-check itself cannot be traced. |
| Oler1961 | OK | P |
| Benavoli2016, Cleghorn2018 | **delete** | Never cited in the supplement. A manual `thebibliography` prints them anyway, and Elsevier requires every listed reference to be cited. |

## 4. Software versions, commits, access dates

- **Main text.** `optA/sw/05_setup.tex` (line 111) states Python 3.11.15, NumPy 2.4.6 and SciPy 1.17.1. These match
  this environment (ENV) and existing upstream tags (`git ls-remote`: cpython v3.11.15, numpy v2.4.6, scipy
  v1.17.1). PyPI upload dates: NumPy 2.4.6 on 2026-05-18, SciPy 1.17.1 on 2026-02-23. Both predate the runs.
- **PyWake.** Version 2.6.20 is consistent across `08_beyond.tex`, the supplement macro `\NFPyWakeVersion`
  (`analysis/mpce_numbers_dir.tex`) and both bibliographies. It is also the installed version.
- **JAX.** Version 0.10.2 appears only in the supplement bibliography. It exists on PyPI, but the script that used it
  is not in the archive.
- **Audit environment.** `SWEVO_rev2/requirements-audit.txt` pins numpy==2.3.5 and pandas==2.2.3. It is labelled
  "Validated standalone audit environment. This is not the missing optimizer environment", so it does not
  contradict the manuscript. A reproducibility package should list the run environment (NumPy 2.4.6, SciPy 1.17.1)
  separately.
- **IEA37 commit.** 267f6e5 is confirmed as HEAD of the IEAWindSystems repository as of 2026-10-04. The docstring of
  `analysis/iea37_model.py` (line 4) names the byuflowlab URL together with commit 267f6e5, but that commit exists
  only in the IEAWindSystems repository; byuflowlab HEAD af88908 has an identical tree. Suggested docstring wording:
  "github.com/IEAWindSystems/iea37-wflo-casestudies, commit 267f6e5 (= byuflowlab af88908 content)".
- **Access dates.** PyWake is accessed Sep. 25, 2026 and IEA37 Sep. 28, 2026. Both are plausible, and IEA37 is
  unchanged as of Oct. 4, 2026. No change is needed.
- **Early-access items.** The only one is Solanki2023, which is still Online First (Section 1). No other entry is
  early access or a preprint-only item, apart from the intentional arXiv citation BartzBeielstein2020.

## 5. Other findings relevant to references and declarations

1. **Reference order: needs action at revision.** The SWEVO guide says references are "numbered in the order they
   appear" (search records of the guide, and Paperpile's style page updated Aug. 2026). In `swevo_back.tex` the
   first 30 entries follow first-citation order, but 45 of 88 do not (from position 31 on: Mosetti1994 is listed 31st
   but first cited 40th, and so on). First-citation order, front to back matter:
   1 Srikakulapu2018, 2 Hou2019, 3 Kusiak2010, 4 Baker2019, 5 Azlan2021, 6 Pookpunt2013, 7 Asaah2021, 8 Mirjalili2017,
   9 Rezk2019, 10 Hooker1995, 11 Derrac2011, 12 BartzBeielstein2020, 13 LaTorre2021, 14 Sorensen2015, 15 Aranha2022,
   16 Wilson2018, 17 Thomas2023, 18 Kennedy1995, 19 Poli2009, 20 Bonyadi2017, 21 Cleghorn2018, 22 Benavoli2016,
   23 Lakens2017, 24 Mladenovic1997, 25 Hansen2001, 26 Solanki2023, 27 Neri2012, 28 Deb2000, 29 Wu2017,
   30 Benavoli2017, 31 Tasgetiren2007, 32 Marinakis2017, 33 Gumaida2019, 34 Feng2015, 35 Castelli2022,
   36 CamachoVillalon2023, 37 Carrasco2020, 38 Guirguis2016, 39 Rodrigues2024, 40 Mosetti1994, 41 Grady2005,
   42 Jensen1983, 43 Katic1986, 44 Bansal2017, 45 Eroglu2012, 46 Eroglu2013, 47 Elkinton2008, 48 Ju2019, 49 Bai2022,
   50 Nagpal2021, 51 Cazzaro2022, 52 IEA37repo, 53 Quaeghebeur2021, 54 Stanley2019, 55 Bastankhah2014, 56 Tao2020,
   57 Quick2023, 58 CriadoRisco2024, 59 Yang2023, 60 Wang2024, 61 Li2025, 62 ParkPark2019, 63 LiZhang2023,
   64 LiRobert2024, 65 Bempedelis2024, 66 Hansen2021, 67 Demsar2006, 68 Holm1979, 69 Campelo2019, 70 Shi1998,
   71 Clerc2002, 72 Trelea2003, 73 VandenBergh2006, 74 Krasnogor2005, 75 Mladenovic2008, 76 Gebraad2016,
   77 Bastankhah2016, 78 Zaharie2002, 79 Deep2007, 80 Storn1997, 81 Kraft1988, 82 Deb1995, 83 DebGoyal1996,
   84 PyWake, 85 Gocmen2016, 86 Carrillo2013, 87 HansenOstermeier2001, 88 Tanabe2014.
   The supplement list is also not in citation order: Eroglu2012 and Eroglu2013 are listed 2nd and 3rd but first
   cited 13th and 14th.
2. **Data availability statement is outdated.** `swevo_back.tex` line 19 and `declarations_content.tex` line 14
   say the archive "does not include the 36,190 original optimization records or all core model, optimizer and
   analysis modules". According to REV3_BRIEF, these files are now in `analysis/`. The statement and supplementary
   Section S-archive need updating; this is not a bibliography change.
3. **Stale source comment.** The comment in `optA/swevo_front.tex` says the abstract is "about 220 words". It is
   242 words (Section 6).

## 6. A7: SWEVO Guide for Authors checklist

The live guide (sciencedirect.com/journal/swarm-and-evolutionary-computation/publish/guide-for-authors) and its
elsevier.com mirror could not be opened (EGRESS_BLOCKED). The requirements below come from search-engine records of
that guide page. Items marked "unverified" are general Elsevier policy or were not found in those records.

| Item | Guide requirement (source) | Manuscript | Result |
|------|---------------------------|------------|--------|
| Abstract | Concise, factual, stand-alone, no references; **max. 250 words** (search record of the guide) | 242 words (LaTeX source tokens), 243 rendered tokens, 246 if the equations are written with spaces as in Word; no citations | **Pass, 4–8 words to spare.** Any addition will exceed the limit. |
| Highlights | Separate editable file with "highlights" in the name; **3–5 bullets, ≤85 characters including spaces each** (search record of the guide) | `SWEVO_highlights.txt`, 5 bullets of 72, 68, 61, 61 and 72 characters | **Pass** |
| Keywords | **1–7 keywords**, English; avoid "and"/"of"; only firmly established abbreviations (search record of the guide) | 6: Benchmarking; equivalence testing; particle swarm optimization; random-sampling control; wake effect; wind farm layout optimization | **Pass** |
| Length | No word or page limit found in the search records | — | **Unverified.** No limit is known; check the live guide. |
| Reference style | Numbers in square brackets, references numbered in order of appearance (search record of the guide; Paperpile style page). Example format "M.M. Stevens, J.H. George, Title, Journal 310 (2005) 1135–1138." Elsevier "Your Paper Your Way" accepts any consistent style at submission. | Numbered `[n]`, consistent IEEE-like style; order **not** by first citation (Section 5.1) | **Action:** reorder (or use `elsarticle-num` with BibTeX). The IEEE-like formatting is acceptable at submission (Your Paper Your Way, from general Elsevier policy; unverified for this journal). |
| Data statement | The journal encourages stating data availability at submission (search record); Research data section in the guide | "Data availability" section present | **Pass on format**; content outdated (Section 5.2) |
| CRediT | CRediT authorship contribution statement expected (seen in published SWEVO articles; guide text not retrieved) | Present. Roles used (Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing – original draft, Validation, Writing – review & editing, Supervision) are all valid CRediT terms. | **Pass** (requirement text unverified) |
| Competing interests | Declaration required for all authors (search record) | "Declaration of competing interest" present, standard Elsevier wording | **Pass** |
| Funding | Elsevier standard | Present ("did not receive any specific grant…") | **Pass** |
| Generative-AI declaration | Section titled "**Declaration of generative AI and AI-assisted technologies in the manuscript preparation process**", placed before the references. Template: "During the preparation of this work the author(s) used [NAME TOOL / SERVICE] in order to [REASON]. After using this tool/service, the author(s) reviewed and edited the content as needed and take(s) full responsibility for the content of the published article." Basic grammar, spelling and reference checkers are exempt. AI used in the research process itself is to be described in the Methods. (Search records of the guide and of Elsevier's generative-AI policy page) | Heading matches exactly and the section is before the references. The wording is not in template form. | **Minor action, suggested below** |

Suggested AI declaration wording. Per A4, it keeps the tools that were actually used. The lead/authors must add any
further tool, including Claude's use in this revision 3:

> During the preparation of this work the authors used Claude (Anthropic) in order to assist with coding and
> manuscript consistency checks, and ChatGPT (OpenAI) in order to assist with the manuscript revision, archive
> audit, document compilation and package preparation. After using these tools, the authors reviewed and edited
> the content as needed and take full responsibility for the content of the published article.

Elsevier's policy also asks that AI used in the research process, such as code generation for the analyses, be
described in the Methods. The authors should decide whether the coding assistance falls under this and, if so, add
one sentence to `05_setup.tex` ("Code, data and environment").
