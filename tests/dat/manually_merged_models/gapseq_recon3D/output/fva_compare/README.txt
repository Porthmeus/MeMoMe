FVA comparison workflow

This folder uses a 3-step pipeline:

1) Run FVA (writes raw FVA CSVs + diagnostics JSON)
   .venv-fva/bin/python \
     tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/01_run_fva_analysis.py \
     --output-dir tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare

   Outputs per model:
   - *_host_fva_baseline.csv
   - *_host_fva_closed.csv (if feasible)
   - *_bacterial_exchange_fva.csv
   - *_bacterial_exchange_nonzero.csv
   - *_diagnostics.json

2) Postprocess (round to 6 decimals, compute deltas/ratios, summary)
   .venv-fva/bin/python \
     tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/02_postprocess_fva_analysis.py \
     --output-dir tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare

   Outputs per model:
   - *_host_fva_comparison.csv
   - *_bacterial_exchange_forced_nonzero.csv
   - summary.txt

3) Plot
   .venv-fva/bin/python \
     tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/03_plot_fva_results.py \
     --output-dir tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare \
     --plot-top-n 30

Note:
- The postprocess step rounds all FVA CSVs to 6 decimals before computing changes.
- If *_diagnostics.json is missing, postprocess fills missing fields with defaults.

4) Microbiome-dependency analysis (loops + dependency lists)
   .venv-fva/bin/python \
     tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare/04_analyze_microbiome_dependency.py \
     --output-dir tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare

   Outputs in: output/fva_compare/microbiome_dependency/
   - loop_only_auto.csv / loop_only_manual.csv / loop_overlap.csv
   - relative_change_scatter_data.csv
   - auto_high_dependency.csv / manual_high_dependency.csv
   - auto_full_dependency.csv / manual_full_dependency.csv
   - top_delta_rel_change.csv
   - summary.txt
