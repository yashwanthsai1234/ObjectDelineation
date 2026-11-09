# Verdict on final_r11_data.md (fresh verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product choice in scope; no catalog lookup, nothing installed.
- engineering:code-review: read `signal.draws` (src/objdel/eval/signal.py:102-107), A1/A4/A-G4 and U11/U51/I22 against the spec text.
- engineering:debug: reproduced the claim chain (pointer -> spec key -> gate input) instead of trusting the review.
- engineering:system-design: traced A1 output keys -> A4 `paired` -> U51 `gate_a` -> A-G4 decision.
- engineering:architecture: checked which file claims to be the rule (A.md header "test writers' contract" vs A4/A-G4).
- engineering:documentation: each verdict cites the source line and the code or measurement.
- huggingface-skills:huggingface-best: no model choice in scope.

## Verdicts
- data-1: CONFIRMED (major). A.md A1 (line 345) says "Output keys: listed in `docs/implementation2.0/specs/A.spec.md` §3.2, plus `gate_ci`" and the A.md header calls A.spec.md "the test writers' contract", while A.spec.md (still "Revision 2") §3.2 says `paired` triples carry "nodule-cluster CIs", `iou_found` "a nodule CI", and §5 item 8 decides Gate A on "its nodule CI lower end", against A4 (A.md:689-693) and A-G4 (A.md:1156) which decide on the patient CI; U11 checks only `[0]==0.0125` and `lo<=d<=hi`, U51 takes one `[d, lo, hi]` triple, I22 checks only presence, so no planned $0 test pins the unit; `signal.draws` today is nodule-only (signal.py:102-107) and the patient CI is 23 % wider on R5 val-select (cluster_ci.json: nodule 0.069, patient 0.085), so the report default can differ; no review log marks it handled; $0 on the Mac, so major not blocking.
- data-2: minor as filed (not promoted): B.spec.md drift is caught by named $0 tests (I07, U49, DG3/DG13 tests, U63); its fix (header precedence rule) rides with data-1.
- data-3: minor as filed: the TOST power bracket is wording; a failed TOST does not stop Gate B1 (B.6.2), no spend changes.
- data-4: minor as filed: carry-over list from r4/r7/r9/r10, not re-argued by the reviewer.

Counts: confirmed 1 (major), refuted 0, downgraded 0.
