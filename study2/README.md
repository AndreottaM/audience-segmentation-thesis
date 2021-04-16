# MM-audiencesegmentation/study2

Contains information on the study investigating the psychological characteristics that underlie audience segments.

## Files

**study2.qsf**: A Qualtrics file that may be imported into their system to reproduce the survey.

The original items used in the mental model scale is shown in **surveyscales-mms.csv** (subscales of items are presented in **surveyscales-mms-with-subscales.csv**). Some minor adjustments to spelling were made (as per Study 1). Some minor adjustments to spelling were made, alongside an additional item to the causal subscale. The original items of the scales were administered in a random order, as Bostrom et al. (2012) does not indicate the original order of administration. The new item was presented as the final item of the causes subscale, for ease of interpretation. Changes are denoted in **surveyscale-changes-mms.R**.

All other survey items are presented in **surveyscales-original.csv**, with references for each scale in **surveyscales-references.docx**/**surveyscales-references.pdf**. Minor changes to make the items relevant for an Australian audience are are noted in **surveyscale-changes.Rmd**. The items presented to participants are shown in **surveysales-new.csv**.

The web application used to collect Q sort data is **qsort-shinyapp.R**. The function **qsconvert.R** was used to transform seven digit codes into seven digit passwords, that are returned to participants upon completion of the Q sort. The statements of the Q sort are shown in **q-statements.csv**. The order of statement presentation is initially randomised using **randomise-statements.R**. The resulting order of presentation is shown in **q-statements-ordered.csv**. The schema used to initialise a SQLite data base is presented in **data/schema.txt**.

The order of presentation of survey scales was counterbalanced using a Latin square design. The counterbalanced orderings are shown in **latin-square.txt**.
