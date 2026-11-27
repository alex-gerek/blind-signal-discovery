# GW170814 cluster-map diagnosis — 27 November 2026

[Full report and maps](../results/gw170814_cluster_maps/20261127T093636Z_698294e8/report.md) · [Predefined protocol](gw170814_cluster_map_protocol_2026-11-27.md) · [Independent CSV/JSON check](../results/gw170814_cluster_maps/20261127T093636Z_698294e8/record_crosscheck.json).

Unchanged H1/L1 event maps reproduce all286qualified raw clusters. Independent graph traversal verifies520components including singletons; all six normalized maps and1,345seed pixels are checked. Original source/input hashes match. No detector modification, calibration or new detection claim.

On the short scale, the strongest reference-associated H1 cluster has20cells and raw score123.225, below192.182; L1 has57cells and464.724, above174.589. Within the actual predefined reference cell region the largest component holds48.79%/76.58% of H1/L1 seed excess. Both sites touch12components (10qualified): simple component count does not establish H1 fragmentation. Medium-scale H1 also has less total seed content. Scattered seeds cannot all be identified as parts of the physical event.

Local temporal normalization differs between sites/regions, but no causal self-normalization loss is established. The report proposes an explicit one-missing-time-frame seed-link variant and a separate paired normalization control. Neither is implemented or accepted. Any detection variant requires full new background normalization/threshold calibration and independent evaluation; regional sums must not be compared with the old cluster threshold as if they recovered H1.

Four-second grid miss results remain unchanged. O3, original strain, other events, detector mathematics, article and CodeQL drafts were not touched. DNA-window idea remains deferred. Research aim: understand the measured mechanism, preserve failures and formulate tests rather than tune a known event into success.
