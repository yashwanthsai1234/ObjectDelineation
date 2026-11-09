
### B1. Final held-out fail rate P(any typed error) at T=600 (mean ± 95% CI, 20 seeds)

| scenario | uniform | diff-filter | learnability | typed-flat | typed-balance | typed-stall | typed-LP |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | 0.085 ± 0.003 | 0.114 ± 0.005 | **0.071 ± 0.002** | 0.073 ± 0.002 | 0.073 ± 0.002 | 0.073 ± 0.002 | 0.073 ± 0.002 |
| S1 default (heterogeneous) | 0.096 ± 0.005 | 0.114 ± 0.006 | **0.070 ± 0.004** | 0.070 ± 0.004 | 0.071 ± 0.004 | 0.070 ± 0.004 | 0.070 ± 0.004 |
| S2 + unlearnable EP (floor 0.3) | 0.111 ± 0.006 | 0.126 ± 0.007 | **0.090 ± 0.004** | 0.093 ± 0.004 | 0.091 ± 0.004 | 0.093 ± 0.004 | 0.093 ± 0.004 |
| S3 + strong forgetting | 0.173 ± 0.007 | 0.216 ± 0.008 | **0.145 ± 0.006** | 0.147 ± 0.007 | 0.147 ± 0.007 | 0.147 ± 0.007 | 0.147 ± 0.007 |
| S4 + strong memorisation | 0.123 ± 0.006 | 0.184 ± 0.007 | **0.116 ± 0.006** | 0.117 ± 0.006 | 0.118 ± 0.006 | 0.117 ± 0.006 | 0.117 ± 0.006 |
| S5 + reward hack | 0.588 ± 0.007 | **0.466 ± 0.013** | 0.594 ± 0.007 | 0.544 ± 0.011 | 0.582 ± 0.008 | 0.529 ± 0.006 | 0.560 ± 0.008 |
| S6 + cross-type transfer 0.5 | 0.067 ± 0.004 | 0.085 ± 0.005 | **0.052 ± 0.003** | 0.053 ± 0.003 | 0.053 ± 0.003 | 0.053 ± 0.003 | 0.053 ± 0.003 |
| S7 small train set (128) | 0.155 ± 0.009 | 0.185 ± 0.012 | **0.154 ± 0.009** | 0.155 ± 0.009 | 0.155 ± 0.009 | 0.155 ± 0.009 | 0.155 ± 0.009 |
| S8 rare type also slow (MH α=0.5, q=0.05) | 0.086 ± 0.004 | 0.102 ± 0.005 | 0.068 ± 0.003 | 0.068 ± 0.003 | 0.068 ± 0.003 | **0.068 ± 0.003** | 0.068 ± 0.003 |
| S9 no diminishing returns (s_cap=∞) | 0.060 ± 0.004 | 0.078 ± 0.006 | **0.038 ± 0.003** | 0.038 ± 0.003 | 0.039 ± 0.003 | 0.039 ± 0.003 | 0.038 ± 0.003 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | 0.240 ± 0.008 | 0.262 ± 0.009 | 0.243 ± 0.008 | 0.258 ± 0.008 | **0.226 ± 0.007** | 0.260 ± 0.008 | 0.253 ± 0.008 |
| S11 slow reward hack (α_h=0.6) | 0.183 ± 0.005 | 0.180 ± 0.008 | 0.166 ± 0.004 | 0.166 ± 0.004 | **0.166 ± 0.004** | 0.166 ± 0.004 | 0.166 ± 0.004 |

### B2. Mean held-out fail rate over training (AUC/T; lower = faster)

| scenario | uniform | diff-filter | learnability | typed-flat | typed-balance | typed-stall | typed-LP |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | 0.170 ± 0.004 | 0.167 ± 0.005 | 0.145 ± 0.004 | 0.136 ± 0.004 | **0.135 ± 0.004** | 0.136 ± 0.003 | 0.136 ± 0.003 |
| S1 default (heterogeneous) | 0.158 ± 0.007 | 0.144 ± 0.006 | 0.121 ± 0.005 | 0.115 ± 0.005 | 0.116 ± 0.005 | 0.115 ± 0.005 | **0.115 ± 0.005** |
| S2 + unlearnable EP (floor 0.3) | 0.171 ± 0.007 | 0.156 ± 0.007 | 0.138 ± 0.005 | 0.134 ± 0.005 | **0.134 ± 0.005** | 0.135 ± 0.006 | 0.134 ± 0.005 |
| S3 + strong forgetting | 0.205 ± 0.008 | 0.203 ± 0.008 | 0.169 ± 0.007 | 0.164 ± 0.007 | 0.165 ± 0.007 | **0.164 ± 0.007** | 0.164 ± 0.007 |
| S4 + strong memorisation | 0.171 ± 0.007 | 0.187 ± 0.007 | 0.150 ± 0.007 | **0.148 ± 0.006** | 0.150 ± 0.006 | 0.148 ± 0.006 | 0.148 ± 0.006 |
| S5 + reward hack | 0.410 ± 0.007 | **0.362 ± 0.009** | 0.412 ± 0.007 | 0.380 ± 0.006 | 0.392 ± 0.007 | 0.379 ± 0.005 | 0.386 ± 0.006 |
| S6 + cross-type transfer 0.5 | 0.121 ± 0.006 | 0.111 ± 0.005 | 0.093 ± 0.004 | 0.089 ± 0.004 | 0.089 ± 0.004 | **0.088 ± 0.004** | 0.089 ± 0.004 |
| S7 small train set (128) | 0.189 ± 0.009 | 0.204 ± 0.012 | 0.178 ± 0.009 | 0.178 ± 0.009 | 0.178 ± 0.009 | **0.177 ± 0.009** | 0.178 ± 0.009 |
| S8 rare type also slow (MH α=0.5, q=0.05) | 0.136 ± 0.005 | 0.125 ± 0.005 | 0.108 ± 0.004 | 0.104 ± 0.004 | 0.105 ± 0.004 | 0.104 ± 0.004 | **0.104 ± 0.004** |
| S9 no diminishing returns (s_cap=∞) | 0.126 ± 0.006 | 0.111 ± 0.006 | 0.087 ± 0.004 | **0.081 ± 0.004** | 0.082 ± 0.004 | 0.081 ± 0.004 | 0.081 ± 0.004 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | 0.294 ± 0.008 | 0.299 ± 0.009 | 0.288 ± 0.009 | 0.283 ± 0.008 | **0.267 ± 0.008** | 0.283 ± 0.008 | 0.280 ± 0.008 |
| S11 slow reward hack (α_h=0.6) | 0.202 ± 0.006 | 0.187 ± 0.007 | 0.177 ± 0.005 | 0.170 ± 0.005 | 0.171 ± 0.005 | 0.170 ± 0.005 | **0.170 ± 0.005** |

### B3. Final held-out worst-type error rate (max over types of per-type error on relevant items)

| scenario | uniform | diff-filter | learnability | typed-flat | typed-balance | typed-stall | typed-LP |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | 0.068 ± 0.003 | 0.101 ± 0.005 | **0.057 ± 0.002** | 0.059 ± 0.002 | 0.059 ± 0.002 | 0.059 ± 0.002 | 0.059 ± 0.003 |
| S1 default (heterogeneous) | 0.337 ± 0.022 | 0.364 ± 0.031 | **0.197 ± 0.016** | 0.200 ± 0.016 | 0.199 ± 0.016 | 0.200 ± 0.016 | 0.200 ± 0.016 |
| S2 + unlearnable EP (floor 0.3) | 0.380 ± 0.012 | 0.388 ± 0.026 | 0.303 ± 0.001 | **0.302 ± 0.000** | 0.302 ± 0.000 | 0.302 ± 0.000 | 0.302 ± 0.000 |
| S3 + strong forgetting | 0.491 ± 0.021 | 0.530 ± 0.028 | **0.359 ± 0.020** | 0.364 ± 0.019 | 0.367 ± 0.020 | 0.363 ± 0.021 | 0.364 ± 0.021 |
| S4 + strong memorisation | 0.393 ± 0.021 | 0.516 ± 0.027 | **0.353 ± 0.021** | 0.361 ± 0.020 | 0.362 ± 0.020 | 0.361 ± 0.020 | 0.360 ± 0.020 |
| S5 + reward hack | 0.553 ± 0.006 | **0.414 ± 0.015** | 0.562 ± 0.007 | 0.511 ± 0.012 | 0.556 ± 0.008 | 0.493 ± 0.005 | 0.529 ± 0.008 |
| S6 + cross-type transfer 0.5 | 0.215 ± 0.017 | 0.253 ± 0.026 | **0.145 ± 0.013** | 0.146 ± 0.013 | 0.147 ± 0.013 | 0.147 ± 0.013 | 0.146 ± 0.013 |
| S7 small train set (128) | 0.473 ± 0.025 | 0.528 ± 0.025 | **0.466 ± 0.025** | 0.471 ± 0.025 | 0.472 ± 0.025 | 0.470 ± 0.025 | 0.471 ± 0.025 |
| S8 rare type also slow (MH α=0.5, q=0.05) | 0.611 ± 0.025 | 0.600 ± 0.028 | 0.451 ± 0.029 | **0.439 ± 0.029** | 0.443 ± 0.028 | 0.441 ± 0.029 | 0.439 ± 0.028 |
| S9 no diminishing returns (s_cap=∞) | 0.255 ± 0.023 | 0.304 ± 0.036 | **0.124 ± 0.013** | 0.126 ± 0.013 | 0.128 ± 0.013 | 0.128 ± 0.013 | 0.127 ± 0.013 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | 0.451 ± 0.000 | 0.463 ± 0.010 | 0.450 ± 0.000 | 0.450 ± 0.000 | 0.450 ± 0.000 | **0.450 ± 0.000** | 0.450 ± 0.000 |
| S11 slow reward hack (α_h=0.6) | 0.336 ± 0.022 | 0.352 ± 0.031 | 0.203 ± 0.017 | **0.201 ± 0.017** | 0.201 ± 0.017 | 0.203 ± 0.017 | 0.202 ± 0.017 |

### B4. Steps to reach the uniform sampler's best held-out fail rate (median over seeds; 'n/m' = seeds that never reached it within 600 steps)

| scenario | target | uniform | diff-filter | learnability | typed-flat | typed-balance | typed-stall | typed-LP |
|---|---|---|---|---|---|---|---|---|
| S0 homogeneous types | 0.085 | 455 (10/20 never) | nan (20/20 never) | 280 | 270 | 260 | 260 | 270 |
| S1 default (heterogeneous) | 0.096 | 480 (10/20 never) | 255 (16/20 never) | 235 | 195 | 195 | 200 | 195 |
| S2 + unlearnable EP (floor 0.3) | 0.111 | 415 (12/20 never) | 235 (16/20 never) | 250 | 240 | 230 | 230 | 235 |
| S3 + strong forgetting | 0.172 | 260 (9/20 never) | 110 (7/20 never) | 130 | 110 | 110 | 110 | 110 |
| S4 + strong memorisation | 0.123 | 375 (10/20 never) | nan (20/20 never) | 210 (5/20 never) | 205 (6/20 never) | 210 (6/20 never) | 205 (6/20 never) | 205 (6/20 never) |
| S5 + reward hack | 0.280 | 70 (9/20 never) | 60 (1/20 never) | 60 (7/20 never) | 60 (1/20 never) | 50 | 60 (1/20 never) | 60 (1/20 never) |
| S6 + cross-type transfer 0.5 | 0.067 | 430 (10/20 never) | 230 (18/20 never) | 210 | 185 | 195 | 185 | 190 |
| S7 small train set (128) | 0.153 | 250 (9/20 never) | 200 (17/20 never) | 180 (7/20 never) | 170 (7/20 never) | 170 (7/20 never) | 165 (8/20 never) | 170 (8/20 never) |
| S8 rare type also slow (MH α=0.5, q=0.05) | 0.086 | 440 (10/20 never) | 240 (17/20 never) | 240 | 205 | 200 | 210 | 205 |
| S9 no diminishing returns (s_cap=∞) | 0.060 | 490 (10/20 never) | 280 (19/20 never) | 215 | 185 | 185 | 190 | 190 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | 0.240 | 460 (9/20 never) | 380 (18/20 never) | 360 (11/20 never) | 220 (15/20 never) | 220 (3/20 never) | 220 (15/20 never) | 200 (13/20 never) |
| S11 slow reward hack (α_h=0.6) | 0.167 | 220 (11/20 never) | 140 (4/20 never) | 155 | 130 | 130 | 130 | 130 |

### B5. Paired differences (same seed = same item population): typed variants minus scalar learnability p(1−p). Negative = typed better; * = 95% CI excludes 0

| scenario | metric | typed-flat | typed-balance | typed-stall | typed-LP | uniform − learnability | diff-filter − learnability |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | final | +0.002 ± 0.000* | +0.002 ± 0.000* | +0.002 ± 0.000* | +0.002 ± 0.000* | +0.014 ± 0.001* | +0.043 ± 0.003* |
| S0 homogeneous types | auc | -0.010 ± 0.000* | -0.010 ± 0.000* | -0.010 ± 0.000* | -0.010 ± 0.000* | +0.025 ± 0.001* | +0.021 ± 0.002* |
| S0 homogeneous types | worst | +0.002 ± 0.001* | +0.001 ± 0.001* | +0.002 ± 0.001* | +0.002 ± 0.001* | +0.011 ± 0.001* | +0.044 ± 0.005* |
| S1 default (heterogeneous) | final | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.026 ± 0.002* | +0.045 ± 0.004* |
| S1 default (heterogeneous) | auc | -0.006 ± 0.000* | -0.005 ± 0.001* | -0.006 ± 0.000* | -0.006 ± 0.000* | +0.037 ± 0.002* | +0.023 ± 0.002* |
| S1 default (heterogeneous) | worst | +0.003 ± 0.002* | +0.003 ± 0.002* | +0.004 ± 0.001* | +0.003 ± 0.002* | +0.141 ± 0.010* | +0.167 ± 0.020* |
| S2 + unlearnable EP (floor 0.3) | final | +0.004 ± 0.001* | +0.002 ± 0.000* | +0.004 ± 0.001* | +0.004 ± 0.001* | +0.021 ± 0.002* | +0.036 ± 0.004* |
| S2 + unlearnable EP (floor 0.3) | auc | -0.004 ± 0.001* | -0.004 ± 0.001* | -0.003 ± 0.001* | -0.004 ± 0.001* | +0.033 ± 0.002* | +0.018 ± 0.002* |
| S2 + unlearnable EP (floor 0.3) | worst | -0.002 ± 0.001* | -0.001 ± 0.000* | -0.002 ± 0.001* | -0.002 ± 0.000* | +0.076 ± 0.011* | +0.085 ± 0.026* |
| S3 + strong forgetting | final | +0.002 ± 0.001* | +0.002 ± 0.001* | +0.002 ± 0.001* | +0.002 ± 0.001* | +0.029 ± 0.001* | +0.072 ± 0.004* |
| S3 + strong forgetting | auc | -0.005 ± 0.001* | -0.004 ± 0.001* | -0.005 ± 0.000* | -0.005 ± 0.000* | +0.036 ± 0.001* | +0.034 ± 0.002* |
| S3 + strong forgetting | worst | +0.004 ± 0.006 | +0.007 ± 0.005* | +0.003 ± 0.007 | +0.005 ± 0.005* | +0.132 ± 0.007* | +0.170 ± 0.015* |
| S4 + strong memorisation | final | +0.001 ± 0.000* | +0.002 ± 0.001* | +0.002 ± 0.000* | +0.001 ± 0.000* | +0.007 ± 0.001* | +0.068 ± 0.003* |
| S4 + strong memorisation | auc | -0.002 ± 0.000* | -0.000 ± 0.000 | -0.002 ± 0.000* | -0.002 ± 0.000* | +0.021 ± 0.001* | +0.037 ± 0.002* |
| S4 + strong memorisation | worst | +0.008 ± 0.002* | +0.009 ± 0.002* | +0.009 ± 0.003* | +0.008 ± 0.002* | +0.040 ± 0.003* | +0.164 ± 0.013* |
| S5 + reward hack | final | -0.050 ± 0.009* | -0.012 ± 0.005* | -0.065 ± 0.006* | -0.035 ± 0.006* | -0.007 ± 0.005* | -0.128 ± 0.012* |
| S5 + reward hack | auc | -0.032 ± 0.004* | -0.020 ± 0.004* | -0.033 ± 0.005* | -0.026 ± 0.004* | -0.002 ± 0.003 | -0.049 ± 0.006* |
| S5 + reward hack | worst | -0.050 ± 0.010* | -0.005 ± 0.006 | -0.069 ± 0.007* | -0.032 ± 0.006* | -0.009 ± 0.005* | -0.148 ± 0.013* |
| S6 + cross-type transfer 0.5 | final | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.015 ± 0.001* | +0.033 ± 0.003* |
| S6 + cross-type transfer 0.5 | auc | -0.004 ± 0.000* | -0.004 ± 0.001* | -0.005 ± 0.000* | -0.004 ± 0.000* | +0.028 ± 0.001* | +0.018 ± 0.002* |
| S6 + cross-type transfer 0.5 | worst | +0.001 ± 0.001 | +0.002 ± 0.001* | +0.002 ± 0.001* | +0.002 ± 0.002 | +0.070 ± 0.006* | +0.108 ± 0.016* |
| S7 small train set (128) | final | +0.001 ± 0.001* | +0.001 ± 0.000* | +0.001 ± 0.001* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.031 ± 0.004* |
| S7 small train set (128) | auc | -0.001 ± 0.000* | +0.000 ± 0.000 | -0.001 ± 0.000* | -0.000 ± 0.000 | +0.011 ± 0.001* | +0.025 ± 0.004* |
| S7 small train set (128) | worst | +0.005 ± 0.002* | +0.006 ± 0.002* | +0.004 ± 0.002* | +0.005 ± 0.002* | +0.008 ± 0.002* | +0.062 ± 0.013* |
| S8 rare type also slow (MH α=0.5, q=0.05) | final | +0.000 ± 0.000 | +0.000 ± 0.001 | -0.000 ± 0.000 | +0.000 ± 0.000 | +0.018 ± 0.001* | +0.034 ± 0.004* |
| S8 rare type also slow (MH α=0.5, q=0.05) | auc | -0.004 ± 0.000* | -0.003 ± 0.001* | -0.004 ± 0.000* | -0.004 ± 0.000* | +0.028 ± 0.002* | +0.017 ± 0.002* |
| S8 rare type also slow (MH α=0.5, q=0.05) | worst | -0.013 ± 0.006* | -0.008 ± 0.006* | -0.010 ± 0.004* | -0.012 ± 0.005* | +0.160 ± 0.011* | +0.149 ± 0.015* |
| S9 no diminishing returns (s_cap=∞) | final | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.001 ± 0.000* | +0.022 ± 0.002* | +0.040 ± 0.004* |
| S9 no diminishing returns (s_cap=∞) | auc | -0.006 ± 0.000* | -0.005 ± 0.000* | -0.005 ± 0.000* | -0.006 ± 0.000* | +0.039 ± 0.002* | +0.024 ± 0.003* |
| S9 no diminishing returns (s_cap=∞) | worst | +0.002 ± 0.003 | +0.004 ± 0.003* | +0.003 ± 0.002* | +0.002 ± 0.003 | +0.131 ± 0.014* | +0.180 ± 0.027* |
| S10 ambiguous type (EP q=0.35, floor 0.45) | final | +0.015 ± 0.002* | -0.016 ± 0.001* | +0.018 ± 0.002* | +0.010 ± 0.002* | -0.003 ± 0.001* | +0.019 ± 0.003* |
| S10 ambiguous type (EP q=0.35, floor 0.45) | auc | -0.005 ± 0.001* | -0.021 ± 0.001* | -0.004 ± 0.001* | -0.007 ± 0.001* | +0.006 ± 0.001* | +0.011 ± 0.002* |
| S10 ambiguous type (EP q=0.35, floor 0.45) | worst | -0.000 ± 0.000* | +0.000 ± 0.000* | -0.000 ± 0.000* | -0.000 ± 0.000* | +0.001 ± 0.000* | +0.013 ± 0.010* |
| S11 slow reward hack (α_h=0.6) | final | -0.000 ± 0.001 | -0.001 ± 0.001 | -0.000 ± 0.001 | -0.000 ± 0.001 | +0.016 ± 0.002* | +0.014 ± 0.004* |
| S11 slow reward hack (α_h=0.6) | auc | -0.007 ± 0.001* | -0.006 ± 0.001* | -0.007 ± 0.001* | -0.007 ± 0.001* | +0.025 ± 0.001* | +0.010 ± 0.002* |
| S11 slow reward hack (α_h=0.6) | worst | -0.002 ± 0.003 | -0.002 ± 0.004 | +0.000 ± 0.004 | -0.001 ± 0.004 | +0.133 ± 0.009* | +0.150 ± 0.020* |

### B6. Overfitting / hacking alarm (train-reward EMA ↑ while some val typed error's 5-eval mean rises above its best earlier 5-eval mean by > max(0.03, 4·SE), for 3 consecutive evals). Action on fire: add the missing guard (EP weight 0.02 → 0.2) and unlearn the hack. Same seeds as B1.

| scenario | sampler | fire rate | median fire step | median step of held-out optimum (no alarm) | final fail w/o alarm | final fail with alarm | best-val-checkpoint (early stop) |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | uniform | 0.00 | nan | 600 | 0.085 ± 0.003 | 0.085 ± 0.003 | 0.087 ± 0.003 |
| S0 homogeneous types | learnability | 0.00 | nan | 555 | 0.071 ± 0.002 | 0.071 ± 0.002 | 0.072 ± 0.002 |
| S0 homogeneous types | typed-balance | 0.00 | nan | 515 | 0.073 ± 0.002 | 0.073 ± 0.002 | 0.074 ± 0.003 |
| S1 default (heterogeneous) | uniform | 0.00 | nan | 600 | 0.096 ± 0.005 | 0.096 ± 0.005 | 0.099 ± 0.006 |
| S1 default (heterogeneous) | learnability | 0.00 | nan | 580 | 0.070 ± 0.004 | 0.070 ± 0.004 | 0.070 ± 0.004 |
| S1 default (heterogeneous) | typed-balance | 0.00 | nan | 545 | 0.071 ± 0.004 | 0.071 ± 0.004 | 0.071 ± 0.004 |
| S2 + unlearnable EP (floor 0.3) | uniform | 0.00 | nan | 600 | 0.111 ± 0.006 | 0.111 ± 0.006 | 0.114 ± 0.006 |
| S2 + unlearnable EP (floor 0.3) | learnability | 0.00 | nan | 545 | 0.090 ± 0.004 | 0.090 ± 0.004 | 0.091 ± 0.005 |
| S2 + unlearnable EP (floor 0.3) | typed-balance | 0.00 | nan | 525 | 0.091 ± 0.004 | 0.091 ± 0.004 | 0.092 ± 0.004 |
| S3 + strong forgetting | uniform | 0.00 | nan | 550 | 0.173 ± 0.007 | 0.173 ± 0.007 | 0.175 ± 0.008 |
| S3 + strong forgetting | learnability | 0.00 | nan | 355 | 0.145 ± 0.006 | 0.145 ± 0.006 | 0.137 ± 0.006 |
| S3 + strong forgetting | typed-balance | 0.05 | 540 | 275 | 0.147 ± 0.007 | 0.147 ± 0.007 | 0.135 ± 0.006 |
| S4 + strong memorisation | uniform | 0.00 | nan | 600 | 0.123 ± 0.006 | 0.123 ± 0.006 | 0.125 ± 0.006 |
| S4 + strong memorisation | learnability | 0.00 | nan | 425 | 0.116 ± 0.006 | 0.116 ± 0.006 | 0.114 ± 0.006 |
| S4 + strong memorisation | typed-balance | 0.00 | nan | 425 | 0.118 ± 0.006 | 0.118 ± 0.006 | 0.116 ± 0.006 |
| S5 + reward hack | uniform | 1.00 | 370 | 100 | 0.588 ± 0.007 | 0.104 ± 0.005 | 0.282 ± 0.007 |
| S5 + reward hack | learnability | 1.00 | 370 | 100 | 0.594 ± 0.007 | 0.096 ± 0.026 | 0.277 ± 0.008 |
| S5 + reward hack | typed-balance | 1.00 | 410 | 100 | 0.582 ± 0.008 | 0.088 ± 0.006 | 0.247 ± 0.007 |
| S6 + cross-type transfer 0.5 | uniform | 0.00 | nan | 600 | 0.067 ± 0.004 | 0.067 ± 0.004 | 0.070 ± 0.005 |
| S6 + cross-type transfer 0.5 | learnability | 0.00 | nan | 545 | 0.052 ± 0.003 | 0.052 ± 0.003 | 0.052 ± 0.003 |
| S6 + cross-type transfer 0.5 | typed-balance | 0.00 | nan | 530 | 0.053 ± 0.003 | 0.053 ± 0.003 | 0.053 ± 0.003 |
| S7 small train set (128) | uniform | 0.00 | nan | 485 | 0.155 ± 0.009 | 0.155 ± 0.009 | 0.154 ± 0.009 |
| S7 small train set (128) | learnability | 0.00 | nan | 360 | 0.154 ± 0.009 | 0.154 ± 0.009 | 0.150 ± 0.009 |
| S7 small train set (128) | typed-balance | 0.00 | nan | 360 | 0.155 ± 0.009 | 0.155 ± 0.009 | 0.152 ± 0.009 |
| S8 rare type also slow (MH α=0.5, q=0.05) | uniform | 0.00 | nan | 600 | 0.086 ± 0.004 | 0.086 ± 0.004 | 0.088 ± 0.004 |
| S8 rare type also slow (MH α=0.5, q=0.05) | learnability | 0.00 | nan | 585 | 0.068 ± 0.003 | 0.068 ± 0.003 | 0.070 ± 0.003 |
| S8 rare type also slow (MH α=0.5, q=0.05) | typed-balance | 0.00 | nan | 540 | 0.068 ± 0.003 | 0.068 ± 0.003 | 0.068 ± 0.003 |
| S9 no diminishing returns (s_cap=∞) | uniform | 0.00 | nan | 600 | 0.060 ± 0.004 | 0.060 ± 0.004 | 0.061 ± 0.004 |
| S9 no diminishing returns (s_cap=∞) | learnability | 0.00 | nan | 575 | 0.038 ± 0.003 | 0.038 ± 0.003 | 0.039 ± 0.003 |
| S9 no diminishing returns (s_cap=∞) | typed-balance | 0.00 | nan | 560 | 0.039 ± 0.003 | 0.039 ± 0.003 | 0.039 ± 0.003 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | uniform | 0.00 | nan | 600 | 0.240 ± 0.008 | 0.240 ± 0.008 | 0.242 ± 0.008 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | learnability | 0.00 | nan | 590 | 0.243 ± 0.008 | 0.243 ± 0.008 | 0.244 ± 0.008 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | typed-balance | 0.00 | nan | 500 | 0.226 ± 0.007 | 0.226 ± 0.007 | 0.227 ± 0.007 |
| S11 slow reward hack (α_h=0.6) | uniform | 0.00 | nan | 360 | 0.183 ± 0.005 | 0.183 ± 0.005 | 0.170 ± 0.005 |
| S11 slow reward hack (α_h=0.6) | learnability | 0.00 | nan | 310 | 0.166 ± 0.004 | 0.166 ± 0.004 | 0.139 ± 0.004 |
| S11 slow reward hack (α_h=0.6) | typed-balance | 0.00 | nan | 280 | 0.166 ± 0.004 | 0.166 ± 0.004 | 0.134 ± 0.004 |

### B7. Generalisation gap (held-out fail − train fail at T) and visit concentration (Gini of per-item visit counts)

| scenario | uniform gap / gini | diff-filter gap / gini | learnability gap / gini | typed-flat gap / gini | typed-balance gap / gini | typed-stall gap / gini | typed-LP gap / gini |
|---|---|---|---|---|---|---|---|
| S0 homogeneous types | 0.024 / 0.09 | 0.023 / 0.77 | 0.025 / 0.26 | 0.027 / 0.30 | 0.027 / 0.30 | 0.027 / 0.30 | 0.027 / 0.30 |
| S1 default (heterogeneous) | 0.025 / 0.09 | 0.023 / 0.80 | 0.030 / 0.33 | 0.032 / 0.37 | 0.033 / 0.37 | 0.032 / 0.37 | 0.032 / 0.37 |
| S2 + unlearnable EP (floor 0.3) | 0.024 / 0.09 | 0.020 / 0.80 | 0.028 / 0.44 | 0.030 / 0.55 | 0.030 / 0.51 | 0.030 / 0.55 | 0.030 / 0.54 |
| S3 + strong forgetting | 0.041 / 0.09 | 0.040 / 0.80 | 0.062 / 0.41 | 0.066 / 0.48 | 0.066 / 0.48 | 0.066 / 0.48 | 0.066 / 0.48 |
| S4 + strong memorisation | 0.088 / 0.09 | 0.068 / 0.82 | 0.092 / 0.22 | 0.093 / 0.23 | 0.094 / 0.23 | 0.093 / 0.23 | 0.093 / 0.23 |
| S5 + reward hack | 0.011 / 0.09 | 0.015 / 0.77 | 0.011 / 0.29 | 0.016 / 0.59 | 0.015 / 0.33 | 0.017 / 0.63 | 0.015 / 0.53 |
| S6 + cross-type transfer 0.5 | 0.018 / 0.09 | 0.017 / 0.81 | 0.021 / 0.28 | 0.022 / 0.31 | 0.022 / 0.31 | 0.022 / 0.32 | 0.022 / 0.31 |
| S7 small train set (128) | 0.135 / 0.04 | 0.110 / 0.49 | 0.137 / 0.12 | 0.138 / 0.13 | 0.138 / 0.13 | 0.138 / 0.13 | 0.138 / 0.13 |
| S8 rare type also slow (MH α=0.5, q=0.05) | 0.022 / 0.09 | 0.021 / 0.81 | 0.031 / 0.33 | 0.033 / 0.37 | 0.033 / 0.37 | 0.033 / 0.37 | 0.033 / 0.37 |
| S9 no diminishing returns (s_cap=∞) | 0.016 / 0.09 | 0.015 / 0.82 | 0.016 / 0.27 | 0.017 / 0.29 | 0.017 / 0.29 | 0.017 / 0.29 | 0.017 / 0.29 |
| S10 ambiguous type (EP q=0.35, floor 0.45) | 0.026 / 0.09 | 0.023 / 0.74 | 0.029 / 0.50 | 0.033 / 0.68 | 0.035 / 0.42 | 0.033 / 0.70 | 0.034 / 0.64 |
| S11 slow reward hack (α_h=0.6) | 0.020 / 0.09 | 0.019 / 0.79 | 0.024 / 0.31 | 0.026 / 0.38 | 0.026 / 0.36 | 0.026 / 0.37 | 0.026 / 0.37 |

### B8. Final held-out per-type error (relevant items), S1 default and S2 unlearnable EP


**S1 default (heterogeneous)**

| sampler | extra_hole | missing_hole | extra_part | missing_part | boundary |
|---|---|---|---|---|---|
| uniform | 0.124 | 0.337 | 0.180 | 0.114 | 0.024 |
| diff-filter | 0.171 | 0.357 | 0.212 | 0.149 | 0.032 |
| learnability | 0.095 | 0.196 | 0.127 | 0.085 | 0.024 |
| typed-flat | 0.095 | 0.199 | 0.130 | 0.087 | 0.024 |
| typed-balance | 0.097 | 0.197 | 0.131 | 0.087 | 0.024 |
| typed-stall | 0.096 | 0.199 | 0.129 | 0.086 | 0.024 |
| typed-LP | 0.095 | 0.199 | 0.127 | 0.087 | 0.024 |

**S2 + unlearnable EP (floor 0.3)**

| sampler | extra_hole | missing_hole | extra_part | missing_part | boundary |
|---|---|---|---|---|---|
| uniform | 0.125 | 0.339 | 0.369 | 0.115 | 0.025 |
| diff-filter | 0.173 | 0.366 | 0.341 | 0.151 | 0.033 |
| learnability | 0.104 | 0.211 | 0.303 | 0.094 | 0.028 |
| typed-flat | 0.111 | 0.218 | 0.302 | 0.101 | 0.030 |
| typed-balance | 0.109 | 0.215 | 0.302 | 0.099 | 0.029 |
| typed-stall | 0.111 | 0.220 | 0.302 | 0.101 | 0.030 |
| typed-LP | 0.112 | 0.219 | 0.302 | 0.102 | 0.030 |

**S8 rare type also slow (MH α=0.5, q=0.05)**

| sampler | extra_hole | missing_hole | extra_part | missing_part | boundary |
|---|---|---|---|---|---|
| uniform | 0.122 | 0.611 | 0.178 | 0.112 | 0.024 |
| diff-filter | 0.169 | 0.600 | 0.212 | 0.149 | 0.031 |
| learnability | 0.094 | 0.451 | 0.126 | 0.085 | 0.024 |
| typed-flat | 0.096 | 0.439 | 0.128 | 0.085 | 0.024 |
| typed-balance | 0.095 | 0.443 | 0.129 | 0.086 | 0.024 |
| typed-stall | 0.094 | 0.441 | 0.128 | 0.085 | 0.024 |
| typed-LP | 0.095 | 0.439 | 0.129 | 0.085 | 0.024 |

**S10 ambiguous type (EP q=0.35, floor 0.45)**

| sampler | extra_hole | missing_hole | extra_part | missing_part | boundary |
|---|---|---|---|---|---|
| uniform | 0.139 | 0.361 | 0.451 | 0.128 | 0.028 |
| diff-filter | 0.189 | 0.420 | 0.450 | 0.171 | 0.043 |
| learnability | 0.145 | 0.323 | 0.450 | 0.131 | 0.039 |
| typed-flat | 0.173 | 0.330 | 0.450 | 0.162 | 0.055 |
| typed-balance | 0.120 | 0.233 | 0.450 | 0.110 | 0.031 |
| typed-stall | 0.176 | 0.342 | 0.450 | 0.162 | 0.057 |
| typed-LP | 0.164 | 0.320 | 0.450 | 0.150 | 0.051 |

### B9. Alarm threshold sensitivity (uniform + learnability samplers, 10 seeds each)

| setting | fire rate, fast hack S5 | median delay after held-out optimum (S5) | fire rate, slow hack S11 | final fail S11 with alarm | fire rate, no-hack scenarios (S0–S4, S6–S10) | worst no-hack scenario |
|---|---|---|---|---|---|---|
| lenient (win 3, z 2, rise 0.02, persist 1) | 1.00 | 135 | 0.90 | 0.101 | 0.850 | S3 (1.00) |
| medium (win 4, z 3, rise 0.02, persist 2) | 1.00 | 230 | 0.15 | 0.164 | 0.040 | S3 (0.20) |
| default (win 5, z 4, rise 0.03, persist 3) | 1.00 | 270 | 0.00 | 0.172 | 0.000 | S0 (0.00) |
| strict (win 6, z 5, rise 0.05, persist 4) | 0.90 | 350 | 0.00 | 0.172 | 0.000 | S0 (0.00) |
