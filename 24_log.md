# Group 24 Analysis Log

2026-09-21 — IC/MP — tested git, created new test repos and played with commit and branching functions. 
2026-09-21 — MP — set up the project repository and report template, added parameters and required libraries, and drafted initial count matrix, sample metadata, and working DESeq2 differential expression analysis code.
2026-09-21 — IC — cloned repo, reviewed assessment structure, adapted Week 6 and 7 count import and metadata code to draft count matrix plus additional key dataframes. Murray's code worked well but we had both been working on the same sections, and after asking unit chair for feedback we decided my code was more aligned to unit content.
2026-09-22 — IC — built the count matrix and sample metadata, checked sample alignment, and merged the completed section into the main report.
2026-09-22 — IC — updated repo with first branch merge to main, included creation of log, proposed git flow structure and commands (to be deleted before zip final analysis).
2026-09-22 — IC — adapted code to filter out unexpressed genes, compute log-CPM for visualisation, and run a principal component analysis to report percentage of variance captured.
2026-09-23 — IC — expanded library-size and filtering checks, centralised the working analysis code, and continued the normalisation and PCA workflow.
2026-09-24 — IC — completed the PCA code and plot, cleaned the section commentary, and updated section 1 plots to order samples by condition.
2026-09-25 — IC — completed significant-gene filtering based on Murray's draft DESeq2 code and plotting and moved low-expression filtering into the normalisation/PCA section to match the assessment workflow.
2026-09-25 — IC/MP — discussed normalisation, CPM visualisation, sequencing-depth correction and the DESeq2 negative-binomial model to clarify interpretation of the normalisation and PCA workflow.
2026-09-28 — MP — added the volcano plot analysis and associated report updates and integrated the significant-gene work into the main analysis.
2026-09-28 — MP — updated PCA commentary and adjusted normalisation and PCA commentary to reflect interpretations of results.
2026-09-28 — IC — completed the draft annotation and GO Biological Process enrichment analysis.
2026-09-28 — IC/MP — reviewed the DESeq2 shrinkage code, identified that a second `lfcShrink()` call was overwriting the apeglm results with normal shrinkage, and removed the redundant assignment.
2026-09-29 — MP — completed the gene-of-interest analysis and interpretation, updated enrichment handling for the zero-enrichment case, and refined interpretation text.
2026-09-29 — MP/IC — investigated the blank up-regulated GO enrichment plot, confirmed that `enrichGO()` returned zero significant pathways, tested the enrichment logic with the non-significant gene set, and reviewed the result together.
2026-09-29 — MP — added summary paragraphs beneath the analysis headers to explain the purpose of each code chunk.
2026-09-29 — MP — added first draft of interpretation for all key sections.
2026-09-29 — IC — reviewed library-size commentary, cleaned comments and reflections in sections 1 and 2, and saved the dispersion plot.
2026-09-30 — MP — refined the library-depth and PCA discussion, supplied the sample_026 interpretation, and created and populated the checked results CSV and helper script.
2026-09-30 — IC/MP — collaboratively reviewed and refined PCA and GO-enrichment interpretation, including the tested-gene universe, annotated-process limitations and biologically cautious wording, replaced some hardcoded values with piped variables.
2026-09-30 — IC — completed final reproducibility updates, including parameterised count paths, saved volcano and heatmap outputs, code and interpretation attributions, and a dynamic volcano caption using the current results.
2026-09-30 — IC — added initials to the setup and session-information chunks and updated checked result values after the final knit with minor adjustments. 
2026-10-01 — MP — corrected the gene-of-interest direction case in the checked results file and helper script.
2026-10-01 — IC — drafted the README and final analysis log for submission and merged the documentation updates.
2026-10-01 — MP - reviewed README and knitted report dress rehearsal.
2026-10-01 — MP — removed unused libraries, reviewed analysis for replicability and accuracy, fixed broken heatmap rendering to save to object, replaced remaning hardcoded values in interpretation with piped variables.
2026-10-01 — MP/IC - final review and sanity checks, knit, ZIP, and submission.
2026-10-02 - MP - replaced boxplot() with ggplot boxplot, tested knitted report on local computer and corrected e- scientific notation to results csv
2026-10-02 - IC - ran final server rmd knit, zipped and checked results from server results, submitted