# Supporting workflows

## Configuration

See [reproduction.md](reproduction.md) for installation, runtimes, licences, native calculations, precision, component experiments, oesophageal inputs, statistics, and diagnostics. Run commands from the repository root:

```bash
python sarc.py [COMMAND] config.json
```

The operation-based commands `analyse`, `evaluate`, `catalogue`, `analyse-edits`, and `integrate` use `operation`, file `inputs`, literal `options`, and a new `output` JSON path:

```json
{"operation":"dawn_freeze","inputs":{"scope":"scope.json","raw_cells":"cells.json"},"options":{"data_class":"PRODUCTION"},"output":"allocations.json"}
```

Input descriptors support `path`, `format` (`json`, `jsonl`, `csv`, `tsv`, `text`, `masked_binary_ledger`), gzip, `paths` lists, `mapping` dictionaries, and JSON-key `select`, except that `catalogue` requires standalone JSON inputs without `select`. Use `masked_binary_ledger` for binary analyses of ledgers also containing continuous outcomes. Other commands use the top-level fields listed below; `statistics` uses its own configuration, not this wrapper. Preserve supplied rosters, seeds, source identities, nulls, unavailable states, and valid-empty outputs.

## Allocation analyses

Command: `analyse`. Question marks mark optional arguments; semicolons separate literal `options`. Exact schemas are linked by the [operation registry](../src/sarc/analysis/run.py).

| Operation | Inputs; options |
| --- | --- |
| `prodigy_freeze` | `scope,carriers` |
| `prodigy_evaluate` | `frozen,labels?` |
| `factorial_decision_freeze` | `protocol,parent` |
| `factorial_decision_evaluate` | `frozen,labels,parent_readout?` |
| `dawn_freeze` | `scope,raw_cells`; `data_class` |
| `dawn_evaluate` | `frozen,selected` |
| `cross_method` | `protocol,pro_freeze,dawn_freeze,pro_labels,dawn_labels` |
| `bridge_freeze` | `protocol,old_freeze,new_freeze,old_domain,new_domain` |
| `bridge_evaluate` | `protocol,frozen,old_labels,new_labels` |
| `native_context_freeze` | `protocol,context` |
| `native_context_merge_labels` | `frozen,documents` |
| `native_context_evaluate` | `protocol,frozen,selected,common_reference` |
| `contextual_order`, `contextual_budgets` | `ranks,eligible,models` |
| `robustness` | `rows,split,bootstrap` |
| `release_consequence` | `rows,split` |
| `source_prior_graph` | `gml,mapping,pathways`; `sample` |
| `source_prior_model` | `scores,statistics,coefficients,eligibility`; `sample,release` |
| `source_prior` | `records` |
| `prior_value_freeze` | `carriers,models,table5_source` |
| `prior_value_evaluate` | `frozen,selected` |
| `prodigy_checkpoint_read` | `requests`; `output_directory,base_directory?,rscript?` |
| `prodigy_retained_hypotheses` | `models,rankings,eligibility,origins,seed_states,results,pathways` |
| `prodigy_hypothesis_census` | `derivative` |
| `prodigy_upstream_selected` | `models,axes,origins,seed_states` |

Freeze eligible rankings before label matching. `MEASURED` calls contain 0/1; `EXPLICIT_NULL` and `ABSENT` contain null. Label sources must match the freeze. PRODIGY identifiers are `root_axis_id,evidence_id,common_domain_id,pathway_dictionary_id,graph_id`; other sources use `source_id`. Decision-adapter names `joint_minus_graph_H0/H1` mean `evidence_G1/G0`.

`source_prior_graph` takes GML as `text` and gene-to-pathway lists. Feed its nested statistics to `source_prior_model` with `{"path":"graph.json","select":["statistics"]}`. `source_prior.records` contains model-operation argument dictionaries. Parsing requires NumPy 1.26.4, pandas 2.2.3, and NetworkX 2.8.8.

The R 4.5.2 [checkpoint reader](../src/sarc/analysis/hypothesis_reader.py) accepts `source_id,mode,path` rows, with mode `checkpoint/result/pathways`. Extract typed `payloads` by source ID for results/pathways. Downstream [hypothesis analysis](../src/sarc/analysis/hypothesis_analysis.py) requires:

- `models`: 85 `sample_ID,model_id,study_donor_id` records.
- `origins`: `source_id,sample,release,checkpoint_source_seed,status,checkpoint`; unavailable checkpoints are null with a reason.
- `seed_states`: complete `sample,release,master_seed,status,mode,checkpoint_source_id,checkpoint_source_seed` records; available hypotheses also require `result_source_id`.
- `rankings`: all 85 x 3 x 5 records with `sample,release,seed,genes,scores`; `eligibility` supplies complete mutation queries.
- `axes`: per-release `genes,samples,normal_samples,sample_origins,source_id`, retaining 34 normal then 85 tumour columns.

Seed modes are `original_full_pipeline`, `deterministic_stage_reuse`, and `original_after_preserved_reuse_assertion`. Reused checkpoints and deterministic aliases are not independent replicates. Census `derivative` is the complete retained-hypothesis output.

## Functional evaluation

Command: `evaluate`; [registry](../src/sarc/evaluation/run.py).

| Operation | Arguments |
| --- | --- |
| `continuous_three`, `continuous_eight` | `selections,endpoints,master_seeds` |
| `matched_return` | `rows,split,rank_rows,eligibility,bootstrap,rank_binding,eligibility_binding` |
| `representation_design` | `records,graphs,selected_samples` |
| `representation` | `observation_document,case_document` |
| `localisation_features` | `baseline,changed,candidates,degs` |
| `localisation_join` | `original,features` |
| `prediction_manifest`, `selected_r2` | `rows` |
| `prediction_metrics` | `rows,predictions`; optional `ids` |

Continuous selections retain model, method, normal pool, seed, policy, status, and fractional weights. Endpoints map sample/model -> gene -> Chronos score or null. Set `options.master_seeds` to `[104729,130363,155921,196613,228017]`; three-model inputs use HCM rather than ACH identifiers. Fixed schemas are in [continuous.py](../src/sarc/evaluation/continuous.py) and [eight_model.py](../src/sarc/evaluation/eight_model.py).

Matched-return inputs require complete ranks, eligibility, paired donor draws, and source bindings. Generate representation designs before evaluating observations. Extract localisation features before joining; edges use `gene1,gene2,combined_score` with integer confidence above 700. Keep applicable `production` defaults enabled.

## Catalogue evaluation

Run `export-catalogue` before `catalogue`. Export configuration contains `operation: "catalogue_cells"`, Boolean `development`, `registry`, `jobs`, and `output`; jobs contain `directory,source_id,cell_ids`. Retain complete score vectors and statuses. [catalogue_export.py](../src/sarc/core/catalogue_export.py) checks 162 development cells across two 36-sample cohorts or 69 protected cells across 396 tumours.

| Catalogue operation | Inputs |
| --- | --- |
| `reference_join` | `primary_rows,secondary_rows,cohort_inputs` |
| `cohort_features` | `samples,old_network,new_network,persistent` |
| `development` | `registry,raw_cells,reference_rows,feature_records` |
| `frozen_choices` | `selection_model,target_features` |
| `protected` | `registry,raw_cells,reference_rows,frozen_choices,frozen_review` |
| `selection_optimism` | `records,populations,biological_units,designs`; optional `ordinary_records` |

Join original NCG6/CancerMine tables using [reference rules](../src/sarc/catalogue/references.py), and build input-only transition [features](../src/sarc/catalogue/features.py). Run development; save its `selection.final_model` object as `selection_model.json` and `review.fitted` object as `frozen_review.json`. Supply these standalone files directly as the corresponding inputs, without `select`. Apply `frozen_choices` to protected input-only features before protected evaluation. For optimism, retain aligned candidate/endpoint records and population maps. [constants.py](../src/sarc/catalogue/constants.py) and [workflow.py](../src/sarc/catalogue/workflow.py) define grids, seeds, and complete support.

## Prediction

Command: `prediction-rounds`. Use Python 3.10.6 and [requirements_prediction.txt](../requirements_prediction.txt). Configuration requires `operation`, `round` (`R1`-`R7`), and a new `output` directory.

| Operation | Additional fields | Output |
| --- | --- | --- |
| `features` (R2/R5/R7) | `slots,molecular`; R5/R7: `pathways` | `features.json` |
| `join` (R2/R5/R7) | `rows,features` | `rows.json` |
| `fit` | `rows`; R6/R7: `partition_order`; R7: `data_id` | `result.json` |
| `summarise` | `rows,result,prior_results` | `summary.json` |

Run features -> join -> fit -> summarise where applicable. R1 uses original rows; R2-R4 the R2 join; R5/R6 the R5 join; R7 its own join from R2 rows. `prior_results` includes preceding rounds, or `{}` for R1. R6/R7 require original 36-sample cohort orders; R7 requires nonempty `data_id`.

Retain all 13,824 development rows and `protected_values_opened: false`. R4 needs exact `event_details.signed_mean_delta`. Feature inputs require complete graph slots and molecular lists; R2 slots include old persistent-baseline `old_rankings`. R5/R7 also need the original triple-tab pathway file. [preparation.py](../src/sarc/prediction_rounds/preparation.py) defines inputs; [adapters.py](../src/sarc/prediction_rounds/adapters.py) fixes grids and partitions.

## Resource comparisons

| Command | Configuration fields |
| --- | --- |
| `compare` | `pairs,output`; optional `comparison` |
| `prepare-temporal` | `sources,baseline_registry,output` |
| `compare-temporal` | prepared `source,output` |
| `canonicalise-string` | compressed-links `source`, Parquet `output` |
| `prepare-string` | `inputs,output` |
| `compare-string` | prepared `source,output` |
| `benchmark` | `baselines,design,output` |
| `benchmark-summary` | `panels,tasks,comparisons,output` |
| `compare-personadrive` | `baseline,followup,samples,settings,output` |

`pairs` contains baseline/follow-up panels with `resource,release,panel_id,members` and optional mappings; `members` are identifier/score pairs. [contracts.py](../src/sarc/resources/contracts.py) defines mapping records.

Temporal preparation precedes comparison: sources map `26.06/26.09` to disease, target, and overall-direct Parquet lists, with the complete [baseline registry](../src/sarc/resources/temporal_prepare.py). For STRING, canonicalise each release, then prepare `baseline_edges,followup_edges,baseline_info,followup_info` and compare.

Benchmarks require panel-to-gene/score baselines and the complete [scenario grid](../src/sarc/resources/benchmark_scenarios.py), retaining distinct `seed/random_seed`. Summarise `task_answers.parquet` and `comparisons.jsonl.gz` with [panel descriptors](../src/sarc/resources/benchmark_summary.py). PersonaDrive uses complete sample-to-gene/score objects, ordered samples, and [comparison settings](../src/sarc/resources/personadrive_comparison.py).

## OpenTargets

```bash
python scripts/opentargets.py download --root /work/opentargets --workers 8
python scripts/opentargets.py analyse --root /work/opentargets --threads 8 --memory-limit 8GB
```

Supply all five Parquet datasets per release under `data/raw/<release>/<dataset>/`: `association_overall_direct`, `association_by_datatype_direct`, `association_by_datasource_direct`, `disease`, and `target`. [support_construction.json](../config/support_construction.json) fixes releases `25.12/26.03/26.06` and source fields; select other configurations with `--config-dir`.

`analyse` runs `harmonise` -> `supports` -> `scores/ranks` -> `sensitivity` -> `robustness`, plus `broader` after supports. Replace `analyse` with one stage only after its dependencies. Retain [release_change_matrix.csv](release_change_matrix.csv), required by sensitivity. Outputs use `data/derived/`; temporary files use `data/interim/`. Use a fresh working directory for a separate analysis.

## Graph interventions

Prepare graphs and controls before analysing ranks. Edge TSVs use canonical `gene1,gene2` and integer `combined_score > 700`; node tables include `preferred_name` and isolates. Use complete cohort mutation matrices and ordered samples.

| Command | Configuration fields |
| --- | --- |
| `prepare-interventions` | `cohorts,transitions,output` |
| `prepare-controls` | `transition,cohort,direction,baseline,followup,persistent_nodes,classes,seeds,output` |
| `assess-residuals` | `inputs.rows,output` |
| `project-residuals` | `inputs.assessment,inputs.scalars,output` |
| `background-contrasts` | `inputs.groups,output` |

Transitions provide `id,baseline,followup,persistent_nodes`. Reverse graphs for retraction; control seeds retain recorded `{seed,random_seed}` pairs 271828/314159/173205. [workflow.py](../src/sarc/interventions/workflow.py) and [analysis.py](../src/sarc/interventions/analysis.py) define graph and residual records.

Command `analyse-edits` supports `edit_slice(document)`, `input_features(baseline,changed,modules,mutations,degs,method)`, and `original_predictors(rows)`. Slice documents contain `cohort,transition,direction,population,graphs,baseline_network_id,references,ranks,bootstrap_indices`. Retain all declared graph slots, complete rankings, five PRODIGY seeds, and 2,000 paired draws. [endpoint_workflow.py](../src/sarc/interventions/endpoint_workflow.py) specifies references and schemas.

## Transfer

Run `select-portfolio`, `evaluate-transfer`, then `transfer-intervals` with separate configurations:

| Command | Configuration fields |
| --- | --- |
| `select-portfolio` | `development,output`; optional `methods` |
| `evaluate-transfer` | `cohort,samples,ranks,support,references,choices,bootstrap,output` |
| `transfer-intervals` | `kind,arguments,output` |

Development rows contain `method,version,score`. Ranks map method -> network -> sample -> `{genes,scores,status}`; preserve complete supports and [reference sets](../src/sarc/transfer/workflow.py). Bootstrap inputs retain `sample_order`, seed 20260926, 2,000 paired indices. Interval `kind` is `selection/transfer/capture`; [intervals.py](../src/sarc/transfer/intervals.py) defines arguments. Measurement keys use `method::version`. Transfer requires both support flags true: `same_ordered_mutation_eligible_sets` and `STRING_values_reprojected_to_fixed_resource_intersection`.

## Integration

Command: `integrate`. Run `describe(catalogue)` before `join_costs(descriptive,non_b_lineage,b_packets,source_ids)`.

The catalogue has `total_rows: 32268`; rows contain `catalog_id,family,original_row_id,source_binding,source_pointer,row`. Supply all 32,088 non-B rows, 379 scenarios, and 180 B policy rows. `source_ids` maps readout/non-B lineage identities; references require nonempty `source_id`. [cost_join.py](../src/sarc/integration/cost_join.py) defines work schemas; [comparison.py](../src/sarc/integration/comparison.py) defines comparison fields. Unknown costs remain unavailable.
