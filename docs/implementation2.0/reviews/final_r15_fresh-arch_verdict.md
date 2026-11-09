# Verdict: final review round 15, fresh-arch

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA product, so no catalog lookup was needed. |
| engineering:code-review | Read the pinned huggingface_hub 2.2.0 download path line by line (`_snapshot_download.py`, `file_download.py`, `_tree_cache.py`, `utils/_xet.py`). This found the per-file HEAD skip that refutes item 1. |
| engineering:debug | Traced one `hf_hub_download` call for a commit-hash revision with a `local_dir` and hf_xet installed: tree cache hit, then `xet_get`, with no `/resolve/` request. |
| engineering:system-design | Re-counted the requests of one rows-root fetch: about 17 paged tree-listing API calls, 1 xet read-token call, and HEAD calls only for the few non-LFS text files. |
| engineering:architecture | Checked whether the one-tar fix improves the design enough to be required. It is a nice-to-have, not a defect fix. |
| engineering:documentation | Each verdict below quotes the source line and the code line. |
| huggingface-skills:huggingface-best | No model choice is in question. Used the HF Hub side only (rate-limit tiers, Xet download path). |

Web check: https://huggingface.co/docs/hub/rate-limits confirms Free user Resolvers 5,000 per 5-minute fixed window,
and "Resolvers" = URLs with a `/resolve/` segment. The numbers in item 1 are right; the request count is not.

## Verdicts

- fresh-arch-1: REFUTED. The item's mechanism ("Each `hf_hub_download` first asks for the file metadata at the file's `/resolve/` URL ... So one rows-root fetch makes at least 16,640 resolver requests") is wrong for the pinned library: in huggingface_hub 2.2.0 [code: experiments/phase2/docker/locks/train.txt:84 `huggingface-hub==2.2.0`, :72 `hf-xet==1.7.0`], `snapshot_download` with a commit-hash revision writes the tree listing to `{local_dir}/.cache/huggingface` [code: .venv huggingface_hub/_snapshot_download.py:251, 282-298], and `_get_metadata_or_catch_error` then "Skip[s] the per-file HEAD call when the file metadata can be rebuilt from a tree listing cached on disk" for every Xet file [code: file_download.py:1736-1746, 1882-1925; local_dir path passes `tree_cache_folder` at :1365-1376], and the Xet download itself uses the once-cached read token, not `/resolve/` [code: file_download.py:2014-2016; utils/_xet.py:161-194]. B1.5 passes exactly a commit hash (`--rows-revision` 40 hex, "`fetch_hf.sh <run> --rows <rows_revision>`", B.md B1.5 "Rows-root download"), and PNG/npz are binary, so the preupload sample (`_commit_api.py:753`) makes them LFS/Xet files. So the fetch costs about 17 tree-listing API calls (the reviewer's own count for 6.1 step 4) plus a handful of HEADs, far under both free-tier buckets; the claimed 10-15 min throttle wait does not occur. What is left is the unpriced transfer minutes, which 00_front 6.2 notes "Why 260 steps" already states ("The rows-root download is not in the measured rows") and already accepts ("If it trips, the adapter is saved and items 3-5 still run"). The `budget_stop_step` label gap is r14 fresh-arch-8 (minor, still open), which this round must not re-raise.
- fresh-arch-2: MINOR (filed as minor; not re-graded). The 00_front 6.2 item 1 row says "rows-root download from HF (pinned revisions)" while its notes say it is not in the measured rows. This is a wording mismatch in one table; it changes no step or value.
- fresh-arch-3: MINOR (filed as minor; not re-graded). The 0.1 min of slack at 27.5 s is a stated choice: 00_front 6.2 notes say what happens when `BudgetStop` trips. Giving the owner a step range is a nice-to-have.

## Counts

confirmed 0, refuted 1 (item 1), downgraded 0. Items 2-3 stay minor as filed.
