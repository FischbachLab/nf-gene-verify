# Gene groups adapted from Santos et al. Figure 1

## Current user-requested override (2026-11-04)

`gene_groups_29.tsv` now follows `data/gene_groups_memo.tsv` for the 22 listed
markers. The other 22 rows preserve their baseline group, color and shape.
Original gene labels are retained: memo `uhpt` matches `uhpT`, and whitespace
around `sdrG` is trimmed. Group-name case is normalized to existing display
headings. These are user-supplied display categories, not new biological
inferences or claims that these categories reproduce the paper exactly.

| Current group | Markers assigned by the memo | Style |
| --- | --- | --- |
| Biofilm | aae, aap, agrA/B/C/D, atlE, bhp, embp, icaA/B/C/D/R, IS256, sarA, sdrF/G/H | Teal circles |
| Virulence | essD, sesI, uhpT | Existing black circles |

`arcA` remains in Virulence because it is not listed in the memo. The ARG
groups and remaining Additional markers retain their earlier assignments and
styles. `data/gene_groups_before_memo.tsv` preserves the earlier configuration;
`gene_groups_29.tsv.memo_update.json` records every applied override.

The latest memo removes `fdh` from its earlier Colonization override. Per the
user's request, `fdh` is restored to **Additional markers**, with the existing
slate-square style. The Colonization panel disappears because it has no
remaining genes. This revision changes no other marker assignment or style.

The earlier paper-derived scheme and its supporting references are retained
below as historical context. The memo overrides the listed assignments.

## Earlier paper-derived scheme

Revised 2026-11-02. This replaces the provisional 29-genome configuration.

The figure uses six antimicrobial-resistance-gene (ARG) color categories,
black circles for virulence/pathogenicity genes and yellow triangles for IS256.
`gene_groups_29.tsv` applies that scheme to the 44 supplied marker labels.
`gene_group_assignment_notes.tsv` records the evidence and qualifications for
every row, distinguishing direct figure matches, name mappings, literature
extensions and markers not classified by the paper.

| Display group | Marker labels | Color and shape |
| --- | --- | --- |
| Beta-lactams | blaZ, mecA | Red squares |
| Quinolones | norA | Purple squares |
| Aminoglycosides | aacA, aadD | Green squares |
| MLSB | ermA, ermB, ermC, mphC, msrA | Orange squares |
| Tetracycline | tetK | Blue squares |
| Other ARGs | dfrC, dfrG, fosB, fusB, fusC | Gray squares |
| Virulence | aae, aap, agrA, agrB, agrC, agrD, arcA, atlE, bhp, embp, icaA, icaB, icaC, icaD, icaR, sarA, sdrF, sdrG, sdrH, sesI | Black circles |
| Additional markers | SCCmec, arcB, arcC, arcD, essD, fdh, uhpT | Slate squares |
| IS256 | IS256 | Yellow triangles |

MLSB abbreviates macrolide, lincosamide and streptogramin B. The figure's gray
"Others" resistance category is called "Other ARGs" here to distinguish it from
the additional markers. Colors reproduce the figure's categorical palette;
these are not claimed to be exact pixel colors from the published image.

## What is directly in Figure 1

The figure explicitly includes blaZ, mecA, norA, aadD, ermA, ermC, msr(A),
mph(C), tetK, fusB, fosB and dfrC in the resistance panel. It also includes
aae, aap, arcA, sdrF, sdrG, sdrH and icaR/icaA/icaD/icaB/icaC in the black-circle
panel. The supplied atlE label is mapped to the paper's atl autolysin label;
this is a name-level mapping rather than a sequence identity check.

The paper's aminoglycoside label AAC(6')-Ie-APH(2'')-Ia is more specific than
the supplied aacA label. aacA is provisionally grouped at the aminoglycoside
family level; the configuration does not assert that it is the same allele
or a complete bifunctional gene. Confirm its reference sequence or annotation
before reporting that exact equivalence.

## Added genes

ermB, dfrG and fusC extend the paper's resistance categories using established
drug classes. agrA/B/C/D, bhp, embp, sarA and sesI extend the black-circle panel
using S. epidermidis studies of regulation, biofilm formation or adhesion.
These additions were not all displayed in Santos Figure 1; they are explicitly
marked as literature extensions in the per-gene notes.

## Markers outside the paper's gene panels

- SCCmec is a cassette-level reported marker, distinct from the mecA gene.
- arcB, arcC and arcD are retained as additional markers; a generic gene label
  does not establish which chromosomal or accessory locus was detected.
- essD is not assigned a S. epidermidis virulence role from its name alone.
- fdh is a formate-dehydrogenase marker, not assigned to an ARG color category.
- uhpT is a transporter-associated marker. Presence of a transporter is not
  equivalent to the resistance-conferring mutations or loss of function studied
  in fosfomycin uptake systems.

arcA remains in the black-circle panel because that placement is visible in
the paper. Its supplied 0/1 call does not establish an ACME-specific arcA copy
or an intact ACME element. No composite ACME call is created from arcA/B/C/D.

All categories describe the reported marker labels. They do not establish an
antibiotic susceptibility phenotype, expression, functional allele, invasiveness,
or biofilm measurement in any of the 29 samples. The supplied calls are unchanged.

## Sources

- S1: [Santos et al. (2025), Figure 1 and results](https://www.nature.com/articles/s41598-025-09061-4).
- S2: [ErmB, CARD ARO:3000375](https://card.csu.mcmaster.ca/ontology/36514), MLSB class.
- S3: [Molecular epidemiology of penicillin-susceptible S. aureus (2022)](https://www.mdpi.com/2076-2607/10/8/1650), dfrG and fusC drug-class annotations. This supports gene-family classification, not a species-specific phenotype prediction.
- S4: [S. epidermidis agr locus polymorphism and pathogenicity (2004)](https://pubmed.ncbi.nlm.nih.gov/15150336/).
- S5: [PIA-, Aap- and Embp-mediated S. epidermidis biofilms (2011)](https://pubmed.ncbi.nlm.nih.gov/21402760/).
- S6: [S. epidermidis orthopedic-device infections (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5625394/), bhp-associated biofilm marker.
- S7: [SarA regulation of S. epidermidis biofilm development (2005)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1065223/).
- S8: [SesI and S. epidermidis invasiveness (2017)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02574/full).
- S9: [S. epidermidis typing and fdh marker (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8106705/).
- S10: [Transporter mutations and fosfomycin resistance in S. aureus (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9525092/). This does not establish a resistance phenotype from uhpT presence in S. epidermidis.
