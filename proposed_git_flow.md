# Possible Git flow for assessment 4

After we've cloned the **24_A4** repo, we start to work from branches that reflect the sections we are working on with roughly standard commit messages.

Commit messages, e.g.,:

- setup project structure
- build count matrix and metadata
- add filtering and PCA
- add DESeq2 analysis
- add shrinkage and significance
- add figures
- add annotation and enrichment
- add interpretation
- final reproducibility check

Branches might be (just an example):

```{text}
main

ic/step1-count-matrix
ic/step2-pca

mp/step3-deseq2
mp/step4-shrinkage

ic/step5-significant-genes
mp/step6-figures
ic/step7-annotation-go
```

## Instructions

Just so I remember the Git commands myself!

First, make sure your local main is current:

```{bash}
cd ~/24_A4

git checkout main
git pull
```

Next, create branch:

```{bash}
git checkout -b ic/step1-count-matrix
```

Then we edit 24_report.Rmd, test that section, and when we reach a checkpoint where we're ready to commit:

```{bash}
git status
```

Then stage the relevant files:

```{bash}
git add 24_report.Rmd 24_log.md
```

Commit:

```{bash}
git commit -m "build count matrix and sample metadata"
```

You might make several commits on the same branch. For example:

```{bash}
git commit -m "add count file sanity checks"
```

then later:

```{bash}
git commit -m "check count matrix against metadata order"
```

So, branch `ic/step1-count-matrix` might contain several commits:

```{text}
build count matrix and sample metadata
add count file sanity checks
check count matrix against metadata order
```

Once Step 1 is working, push the branch:

```{bash}
git push -u origin ic/step1-count-matrix
```

Then on GitHub we can create a pull request into main, ideally have eachother look over it, and merge it.

After it's merged:

```{bash}
git checkout main
git pull
```

Then create the next branch:

```{bash}
git checkout -b ic/step2-pca
```

We can also log messages such as "2026-09-21 — IC — cloned repo, reviewed assessment structure, adapted Week 7 count import and metadata code" in the log file, to complement the commits.

This is what my ~ folder looks like:

/home/[student_id]/
├── 24_A4/
│   ├── 24_report.Rmd
│   ├── 24_report.html
│   ├── 24_results.csv
│   ├── README.md
│   ├── 24_log.md
│   ├── figures/
        ├── PCA_plot.png
        ├── dispersion_plot.png
        ├── volcano_plot.png
        └── heatmap.png
└── 24_A4_working/
    ├── outputs/
    └── scripts/


