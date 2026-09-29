# SWEVO (Swarm and Evolutionary Computation, Elsevier, ISSN 2210-6502): submission checklist

Researched 2026-09-29. **Caveat:** sciencedirect.com, elsevier.com and every mirror were blocked by the network egress proxy, so no Elsevier page could be opened directly. Everything below comes from web-search result snippets of the official Guide for Authors (GfA) and Elsevier policy pages, plus standard Elsevier practice. Labels used below:
- **[V]** means verified from a snippet of the SWEVO GfA or Elsevier policy page.
- **[E]** means the item is Elsevier-wide standard and very likely applies, but I could not confirm it in the SWEVO text.
- **[A]** means assumed or unknown. Check it in Editorial Manager before submitting.

Official GfA: https://www.sciencedirect.com/journal/swarm-and-evolutionary-computation/publish/guide-for-authors

## Article type and review
- [ ] **Full-length research article.** The journal publishes original research articles, surveys of timely topics, and novel applications [V]. Choose "Research paper" / "Full Length Article" in Editorial Manager. The exact menu label is [A].
- [ ] **Review model is single-anonymized** [V, snippet]: reviewers know who the authors are. The editor screens the paper first, then sends it to at least 2 reviewers [V].
  - This means **no anonymized manuscript is needed**. Authors, affiliations and acknowledgements can stay in the main file.
  - A separate title page is not required for anonymity. EM may still ask for the title-page data in metadata fields [E].

## Front matter
- [ ] **Abstract:** about 250 words [V, weak]. A snippet says "short abstract of about 250 words". Keep it to 250 words or fewer, with no references and no undefined abbreviations [E].
- [ ] **Keywords:** 1–7 [V]. Avoid multi-word keywords joined by "and" or "of" [V].
- [ ] **Highlights:** 3–5 bullets, each 85 characters or fewer including spaces [V].
  - Upload them as a **separate editable file** named with "Highlights" in the file name [E].
  - Recommended at submission [E]. Treat them as required.
- [ ] **Graphical abstract:** optional but encouraged [V].
  - If provided, upload it as a separate file [E].
  - Elsevier's guidance: at least 531×1328 px (h×w), readable at 5×13 cm, TIFF/EPS/PDF/MS Office format [E].

## Format and LaTeX
- [ ] **Template:** LaTeX submissions are accepted and the journal's LaTeX template is encouraged [V]. That template is **elsarticle** (CTAN / Overleaf "Elsevier article") [E].
  - Suggested setup for review: `\documentclass[preprint,review,12pt]{elsarticle}`. `review` gives double spacing. Use `3p`/`5p` for a journal-like preview [E].
  - Upload the .tex source, the .bib/.bbl and the figures, plus a compiled PDF [E/V: ".tex extension"].
- [ ] **Columns:** Word files must be single-column. Double-column layout is allowed only for LaTeX submissions [V].
- [ ] **"Your Paper Your Way" at first submission** [E, weak]: format-free submission is likely accepted, i.e. any consistent reference style and one PDF. The SWEVO page's YPYW wording was not seen directly. Journal-style formatting is needed at revision/acceptance.
- [ ] **Page or word limits:** none found [A/unknown]. Snippets show no explicit limit. Keep the paper concise, roughly 25–35 pages in single-column preprint format.

## References
- [ ] **Numbered style** with square brackets [1], [2] [V, snippet plus CSL styles].
  - Use `\bibliographystyle{elsarticle-num}` (or `elsarticle-num-names`) with `\usepackage[sort&compress]{natbib}` or `\biboptions{sort&compress}` [E].
- [ ] Every in-text citation must be in the reference list and vice versa. References cited in the abstract must be given in full [V].
- [ ] Include DOIs and cite datasets and software in the reference list [E].

## Figures, tables and supplements
- [ ] **Figures:** use vector formats (EPS/PDF) for line art. Raster images need ≥300 dpi for halftones, ≥500 dpi for combination art and ≥1000 dpi for line art [E].
  - Figures must be readable in both colour and greyscale. Colour online is free [E].
  - Each figure needs a caption. Number figures in the order they are cited [E].
- [ ] **Tables:** editable text, not images. Number them consecutively, and put the caption above the table and notes below it [E]. Avoid vertical rules [E].
- [ ] **Supplementary material:** submit it as a separate file (PDF, zip, code) with a short caption. It is published as supplied, without typesetting [E].
- [ ] **Code/data:** Mendeley Data is supported for SWEVO [V: journal Mendeley Data page exists]. You can deposit code and algorithms there.

## Required statements (end of manuscript, before References)
- [ ] **Data availability statement.** SWEVO reportedly follows Elsevier research-data **Option C** [V, weak]: deposit the data, cite and link it, or state why it cannot be shared. EM also has a data-statement step [E].
  - Suggested wording: "Data and code are available at [URL/DOI]." or "Data will be made available on request."
- [ ] **Declaration of competing interest** [V].
  - Generate it with Elsevier's Declaration of Interests tool, pick "I have nothing to declare" if that applies, and upload the resulting .doc/.docx at the attach-files step [V].
  - Repeat the statement in the manuscript [E].
- [ ] **Funding statement** [V]: list the funding sources and describe each sponsor's role. If there is no funding, write: "This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors." [E, standard Elsevier wording]
- [ ] **CRediT author statement** [E]: roles such as Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Writing – original draft, Writing – review & editing, Visualization, Supervision. Not verified as mandatory for SWEVO. Include it anyway because it is expected at revision.
- [ ] **Declaration of generative AI** [V, Elsevier-wide]. Add it only if AI was used for writing; basic grammar and spell checks are exempt. Current heading (2025 Elsevier update):
  - **"Declaration of generative AI and AI-assisted technologies in the manuscript preparation process"**. The older heading was "...in the writing process".
  - Place it at the end, before References.
  - Wording: *"During the preparation of this work the author(s) used [NAME OF TOOL / SERVICE] in order to [REASON]. After using this tool/service, the author(s) reviewed and edited the content as needed and take(s) full responsibility for the content of the published article."*
  - AI tools cannot be listed as authors.
  - AI-generated or altered images must also be disclosed in the caption [V, Elsevier 2025 policy].

## Submission logistics
- [ ] **Cover letter** [E/A]: EM has a cover-letter field. State novelty, fit with SWEVO's scope, and that the paper is not under review elsewhere. If this is a resubmission or revision, attach the response-to-reviewers letter.
- [ ] **Suggested reviewers** [A]: EM usually asks for 3–5 potential reviewers with institutional emails. Avoid conflicts of interest. Whether SWEVO requires them is unverified.
- [ ] **ORCID** [E]: recommended for all authors and linked in EM. It is typically required for the corresponding author.
- [ ] **Preprint policy** [E]: Elsevier allows preprints (e.g. arXiv) at any time. Disclose the preprint at submission. Once published, add a link to the DOI.
- [ ] **Open access / APC** [V, snippet]: hybrid journal (subscription or gold OA). The APC is about **USD 3,160** excluding taxes. Recheck on the journal page; institutional agreements may cover it.

## Editors-in-Chief
- **Swagatam Das** (Indian Statistical Institute, Kolkata) and **P. N. Suganthan**, founding co-EiCs [V from editorial-board snippets]. Suganthan is now at the KINDI Centre, Qatar University [V, ORCID/bio]; the board page may still show NTU Singapore.
- Their current EiC status as of 2026 is **not directly verified** because the board page was blocked. Confirm on https://www.sciencedirect.com/journal/swarm-and-evolutionary-computation/about/editorial-board before addressing the cover letter. A safe salutation is "Dear Editors-in-Chief".

## Sources (search-result snippets only)
- SWEVO GfA: https://www.sciencedirect.com/journal/swarm-and-evolutionary-computation/publish/guide-for-authors
- SWEVO GfA mirror PDF: https://apps.lib.whu.edu.cn/ensci/editnew/upfile/2210-6502.pdf
- SWEVO editorial board: https://www.sciencedirect.com/journal/swarm-and-evolutionary-computation/about/editorial-board
- Elsevier GenAI policy for journals: https://www.elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals
- Elsevier updated GenAI policy announcement: https://www.elsevier.com/connect/updated-generative-ai-policies-for-journals-supporting-responsible-use-while-protecting-trust
- SWEVO on Mendeley Data: https://data.mendeley.com/journal/22106502
- SWEVO citation style: https://paperpile.com/s/swarm-and-evolutionary-computation-citation-style/
