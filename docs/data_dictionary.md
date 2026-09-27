# Data dictionary

| Source | Table | Relevant fields | Availability |
| --- | --- | --- | --- |
| RunParameters XML and SampleSheet | `run_lots` | RunID, RunName, Panel, RunDate, Instrument, software and consumable identifiers | Instrument run folder processing only |
| InterOp sample index summary | `percent_index` | Sample_ID, Project, Index_1_I7, Index_2_I5, Percent_Reads_Identified | Instrument run folder processing only |
| InterOp run summary | `summary` | Density, PF_Percent, Percent_Q30_R1, Percent_Error_Rate_R1, Percent_Aligned | Instrument run folder processing; archived output also supplied as CSV |

The archived CSV has one row per run. `Density` is in K/mm². Percent fields retain their SAV export values. The supplied export labels the NextSeq 550 family as `NextSeq`; derived platform files use `NextSeq 550` to make the family clear. `Percent_Aligned` corresponds to the dashboard trace labeled “% PhiX Aligned”; the actual alignment target depends on local run configuration. Empty measurements remain empty.

The supplied CSV contains operational RunID and RunName values. It has no patient fields.
