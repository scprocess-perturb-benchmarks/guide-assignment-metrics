# guide-assignment-metrics

Evaluation module for the [guide-assignment](https://github.com/tmp-perturb/guide-assignment)
Omnibenchmark. Active public entrypoints and outputs follow the benchmark
methodology; archived compatibility metrics remain callable but are not wired in
the plan.

| Entrypoint | Role |
|---|---|
| `run.py` (default) | per-lineage main metrics: `reference_agreement`, `knockdown`, `differential_expression_discovery`, difficulty-stratified counterparts, and `method_resources`; `construct_set_agreement` is supplementary |
| `collectors_run.py` | supplementary `pairwise_top1_agreement` across assignment configurations |
| `aggregate_run.py` | writes separate `main_metrics` and `supplementary_metrics` tidy tables (TSV, plus Parquet when available) |

- Env: `assignment_metrics`
