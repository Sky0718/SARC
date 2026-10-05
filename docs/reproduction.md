# Reproduction guide

## Setup

Run `python sarc.py COMMAND CONFIG.json` from the repository root. Paths are configuration-relative; executable names may resolve on `PATH`. Use new output directories/files. Data and external methods are supplied separately; update configuration paths to their corresponding complete inputs.

Native R calculations require R 4.5.2, Matrix 1.7.4, `plyr`, `data.table`, `jsonlite` and upstream dependencies. PRODIGY additionally uses DESeq2 1.50.2, igraph 2.3.3, PCSF 0.99.1, MASS 7.3.65 and mixtools 2.0.0.1. Supply `rscript`, ordered `library_paths` and a package `versions` map. Oesophageal preparation requires edgeR 4.8.2. Molecular preparation/PersonaDrive use NumPy 1.26.4 and pandas 2.2.3; PersonaDrive also requires NetworkX 2.8.8 and tqdm 4.66.5. `python_executable` selects that environment. Precision checking uses SciPy; statistics require Python 3.10+ and NumPy 1.22+. Resource dependencies are in [requirements.txt](../requirements.txt). Keep compatible numerical environments within comparisons.

## Inputs

Retain complete populations: colorectal 85 models/83 donors; oesophageal 59 models/58 donors. Models are averaged within donor, then donors equally. Primary seeds are `[104729,130363,155921,196613,228017]`; deterministic methods retain aliases, not independent replicates.

CSV matrices have unique ordered gene rows and sample columns, with literal gene IDs first. Mutations are binary. Network TSV columns are `gene1,gene2,combined_score`, containing integer `700 < score <= 1000`, unique undirected non-self edges and preserved row order.

Allocation JSON contains `models` with `model_id,sample_id,donor_id,query` and `allocations` rows with `model_id,cell,budget,weights`. `state` may replace `cell`; four cells are `G0H0,G1H0,G0H1,G1H1`, with underscores also accepted. Weights are gene-to-rational mappings or unique gene/value pairs. PRODIGY rows require `seed` and the complete emitted `menu`; all five seeds must be present. `SUCCESS_EMPTY` requires an empty menu. Capacity is the emitted-menu/query intersection capped by budget; do not refill short menus. Deterministic aliases must agree. Degree uses graph-wide degree restricted to the full query, including eligible zero-degree genes. Unavailable predictions are not zeros.

Label JSON has `calls` with `sample_id,gene,state,value`. `MEASURED` requires integer 0/1. `ABSENT,EXPLICIT_NULL,IDENTITY_CONFLICT,CONFLICTING_CALLS` require null. All allocated identities require explicit calls and repeated identities must agree. Diagnostics additionally require `source_binding.source_id`. Keep candidate eligibility independent of label availability.

Colorectal category JSON contains `sets.organoid` (751 symbols), `sets.cell_line` (1,121) and `roles[gene].cosmic_tokens`. Oesophageal scope supplies `categories.core` and `categories.role` for every query gene. Core and gene-role classes form separate partitions; preserve literal source membership.

## Core calculations

```bash
python sarc.py prepare-molecular molecular.json
python sarc.py prepare common_domain.json
python sarc.py compute native_job.json
python sarc.py export-native export.json
```

`prepare-molecular` requires `contract,source_root,observed_genes,mutation,mutation_sample_map,network,expected_models,output`; the source-column/transform contract follows [molecular.py](../src/sarc/core/molecular.py). `prepare` takes two ordered paths each for `counts,mutations`, plus `models,expected_models,expected_genes,output`; optional inherited evidence preserves source/model/seed identities. See [prepare.py](../src/sarc/core/prepare.py).

Every native job declares `method,mode,models,expected_models,context,output,seeds`. Method-specific fields:

- PRODIGY: `original_sources,counts,mutations,pathways,sample_origins,alpha,cells,tasks`. Origins use literal `normal`/`tumor`. Tasks contain `cell,sample,master_seed,patient_seed,network_id,network`, covering every model/cell/seed. Preserve recorded patient seeds. Factorial order is `G0_H0,G1_H0,G0_H1,G1_H1`; tasks supply inherited `explicit_evidence`, without refitting it. Its available object uses `schema="cdf_explicit_evidence_v1"`, ordered genes and little-endian `binary64_le`; see [prodigy.R](../src/sarc/core/r/prodigy.R). Primary `alpha=0.05`.
- DawnRank native: `tumour,normal,mutation,network,network_id,normalisation_names,normalise_source,expected_normals`. Name CSV columns are `role,canonical,original`. Factorial inputs instead require `common_genes,mutation,edges_G0,edges_G1,inherited_H0,inherited_H1,native_genes_H0,native_genes_H1,expected_genes`. Set `mu=3` and explicit `epsilon,max_iterations`.
- PersonaDrive, native-only: `loader_source,ranking_source,pathways,mutations,outliers,network,network_id,cancer,dataset,pipeline_sample_order,pipeline_to_canonical`. Mutations are gene-by-sample; outliers sample-by-gene. Preserve full-cohort construction and bijective pipeline/canonical identities. See [persona.py](../src/sarc/core/persona.py).

Export PRODIGY with `operation="prodigy_carriers"` and `directory,scope,cell_ids,source_id,output`; `cell_ids` rows contain `model_id,graph,evidence,signature_id`. Keep complete emitted ranks. Other exports are defined in [export.py](../src/sarc/core/export.py) and [selections.py](../src/sarc/core/selections.py).

The main chain is molecular preparation -> all-five-seed PRODIGY/export -> strict DawnRank below -> `diagnostics_prodigy.json` -> `diagnostics_normalise.json` -> paired updates, donor draws, endpoint envelopes, categories and category sensitivity. Supporting analyses are listed in [analyses.md](analyses.md).

## Numerical precision

```bash
python sarc.py precision config/precision_first.json
python sarc.py verify-precision config/precision_verify.json
```

[precision_second.json](../config/precision_second.json) specifies the second transition. Prepared RDS requires ordered `genes`, `models` as sample names or a data frame with `sample_id`, named binary `mutation`, `graphs$G0/G1$transition` and `$damping`, and `evidence$H0/H1` matrices or lists containing `values,evaluable`. Supply actual graph `degree` or `adjacency` for degree export. The roster JSON contains `models` with all three identity fields. [batch.R](../src/sarc/precision/batch.R) validates axes.

Current settings are `epsilon=1e-12,max_iterations=10000`; the original-threshold comparison uses `1e-4,100`, intermediate `1e-6`/`1e-8` settings use 10,000. Stopping and cutoff certification remain separate. Outputs include `allocations.json`, `degree_allocations.json` when available, `cutoff_certificates.csv` and complete binary64 arrays. Arrays are little-endian, column-major, with `genes.txt`/`models.txt` axes. [solver.R](../src/sarc/precision/solver.R) and [certificate.R](../src/sarc/precision/certificate.R) supply the shared implementation.

## Component comparisons

```bash
python sarc.py precision config/precision_components.json
python sarc.py mechanism config/mechanism_freeze.json
python sarc.py mechanism config/mechanism_evaluate.json
python sarc.py mechanism config/mechanism_score.json
```

The eight states `Pa_Db_Hh` combine propagation, damping/restart and evidence on the first colorectal domain. `reuse_directory` supplies matching strict diagonal states. Freeze complete allocations/coefficient requests before selected-label linkage; evaluation requires exactly those calls and fixed core sets. Outputs retain all 22 contrasts, unknown states and donor weights. The score-reference step consumes component scores/categories, preserves negative reference values and exact cutoff ties, and outputs signed allocation terms, query-gene margins and patterns. See [workflow.py](../src/sarc/mechanism/workflow.py) and [score_reference.py](../src/sarc/mechanism/score_reference.py).

## Oesophageal comparison

[oesophageal_sources.json](../config/oesophageal_sources.json) records all models, 403 GTEx gastro-oesophageal-junction references, 118 VCFs, source URLs and files. Populate those paths without changing rosters or exclusions. Retain tumour HTSeq counts, normal gene-read counts and complete current identifier metadata. Preparation names both native evidence domains, graphs and original `DawnNormalize.R`.

```bash
python sarc.py oesophageal config/oesophageal_inputs.json
python sarc.py oesophageal config/oesophageal_prepare.json
python sarc.py precision config/oesophageal_compute.json
python sarc.py oesophageal config/oesophageal_allocations.json
python sarc.py oesophageal config/oesophageal_scope.json
python sarc.py oesophageal config/oesophageal_freeze.json
python sarc.py oesophageal config/oesophageal_evaluate.json
```

Preparation normalises all 462 columns together, applies native-domain normalisation, then projects to 7,388 common genes. Allocation requires certified cutoffs. Outputs are `results/esca_allocations/allocations.json`, `results/esca_scope/scope.json` and `results/esca_evaluation`, including labels, complete ledgers, `composition.json` and `DONOR_DRAWS.npy`.

## Statistics and diagnostics

Run each configuration with `python sarc.py statistics config/NAME.json`:

| Configuration stem | Output |
| --- | --- |
| `statistics_paired_crc` | Absolute/direct degree-relative updates, seed opposition |
| `statistics_draws_crc` | Donor order and 10,000 paired draws |
| `statistics_endpoint_crc` | Eight-coordinate update family |
| `statistics_categories_crc`, `statistics_categories_esca` | Complete category ledgers/contributions |
| `statistics_category_sensitivity_crc`, `statistics_category_sensitivity_esca` | Quantity-fraction families |
| `statistics_signed_accounting` | Signed donor masses |

Statistics use `results/crc_precision/degree_allocations.json`, not menu-restricted degree output. Colorectal draws use `donor_order.json` and `draw_indices.json`; oesophageal evaluation supplies `composition.json` donor IDs and `DONOR_DRAWS.npy`. Reuse each context's complete plan: 10,000 integer-index rows of length 83 or 58, aligned to lexicographic donor order.

Category sensitivity requires fully measured focal allocations. Ratios use pooled positive/returned masses and are recomputed per draw. Zero-allocation donors remain included; undefined denominators are recorded, never dropped. Update families have eight colorectal or four oesophageal coordinates; category families have 12 or eight, respectively, and receive separate calibration. Exact constants retain point envelopes; nonconstant zero-scale coordinates make the planned standardised family unavailable. Unknown-label bounds and donor-composition envelopes are different outputs.

Run `python sarc.py diagnostics config/NAME.json` for:

| Configuration stem | Purpose |
| --- | --- |
| `diagnostics_prodigy`, `diagnostics_normalise` | Native export and combined model ledger |
| `diagnostics_rank`, `diagnostics_influence` | Budgets, seeds, turnover, label/gene omissions |
| `diagnostics_native_menus`, `diagnostics_degree_menu` | Emitted-menu/full-query degree comparisons |
| `diagnostics_menu_records`, `diagnostics_opportunity` | Menu-conditioned action overlap |
| `diagnostics_offered_eight` | Offered-budget continuous totals |

Preserve all budgets 1/5/10/20 and complete menus; `normalise`/`rank` require the complete four-budget matrices and certified DawnRank ten-slot allocations. Label support must include degree and every requested budget, not only native top-ten candidates. Degree inputs require ordered genes and graph TSVs. Full schemas are in [statistics/workflow.py](../src/sarc/statistics/workflow.py) and [diagnostics/workflow.py](../src/sarc/diagnostics/workflow.py).

## Licences

Project code uses Apache License 2.0. Obtain external methods separately and retain their notices:

- [TARGET-SL](https://github.com/RhysGillman/TARGET-SL/tree/ec057f299b3b7017e743660b3a96bcaab9e7fb4f), revision `ec057f299b3b7017e743660b3a96bcaab9e7fb4f`, GPL-3.0. PRODIGY requires `data/PRODIGY_pathways.rds` and these `scripts/PRODIGY/` sources in order: `PRODIGY.R`, `run_single_PCSF.R`, `get_pathway_network.R`, `get_enrichment_pvalue.R`, `get_enriched_pathways.R`, `get_diff_expressed_genes.R`, `get_DEGs.R`, `analyze_PRODIGY_results.R`. DawnRank requires `scripts/DawnRank/DawnNormalize.R`; [DawnRank 1.2 lineage](https://github.com/MartinFXP/DawnRank/tree/e34fefd2e0d52ff6fa8b3384cbf6d24eac4766fd) declares GPL.
- [PCSF 0.99.1](https://github.com/IOR-Bioinformatics/PCSF/tree/4e5f2707fdd096b760e5bb4ca60938a3c9836d78), revision `4e5f2707fdd096b760e5bb4ca60938a3c9836d78`, including solver and licence.
- [PersonaDrive](https://github.com/abu-compbio/PersonaDrive/tree/880de35845180db27a8edf068c5a6bd5c95f333a), revision `880de35845180db27a8edf068c5a6bd5c95f333a`: `constructing_PBNs.py`, `PersonaDrive.py`, `data/kegg_pathways_v1.txt`. No licence file is supplied at that revision; obtain permission separately.

Retain [RBO attribution](../licenses/RBO.txt). Data access does not itself grant redistribution rights; obtain specified releases through official distributions.
