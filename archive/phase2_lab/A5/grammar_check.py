"""A5 sandbox: does the D10 regex compile in xgrammar and llguidance for our three tokenizers, accept every
SFT target (canonical tokenization + the model's stop token), reject malformed texts, and what does it cost?

Usage: venv/bin/python -I scripts/grammar_check.py data/targets_texts.json [--n N] [--models a,b]
"""
import json
import random
import re
import statistics as st
import sys
import time

import llguidance as llg
import llguidance.hf
import llguidance.numpy as lln
import numpy as np
import xgrammar as xgr
from transformers import AutoTokenizer

MODELS = {
    "smolvlm-256m": ("HuggingFaceTB/SmolVLM-256M-Instruct", [49279]),  # generation_config eos <end_of_utterance>
    "lfm2.5-vl-450m": ("LiquidAI/LFM2.5-VL-450M", None),  # tokenizer eos (<|im_end|>)
    "qwen3-vl-2b": ("Qwen/Qwen3-VL-2B-Instruct", [151645, 151643]),
}


def malformed(t: str) -> dict[str, str]:
    """Texts the strict parser rejects, built from a valid target."""
    i = t.index("[[") + 2
    j = t.index(",", i)
    first = t[i:j]
    return {
        "leading_zero": t[:i] + "0" + first + t[j:],
        "range_1001": t[:i] + "1001" + t[j:],
        "space": t.replace(",", ", ", 1),
        "no_holes": t.replace(',"holes":[]', "", 1),
        "keys_order": t.replace('"outer"', '"OUTER"', 1),
        "trailing_text": t + " ok",
        "float": t[:i] + first + ".5" + t[j:],
        "negative": t[:i] + "-" + first + t[j:],
        "two_points": re.sub(r'"outer":\[\[.*?\]\],"holes"', '"outer":[[1,2],[3,4]],"holes"', t, count=1),
        "truncated": t[: len(t) // 2],
        "unicode_digit": t[:i] + "١" + t[j:],
    }


def xgr_run(compiled, tok_info, ids, stop, bitmask, timings=None):
    m = xgr.GrammarMatcher(compiled, override_stop_tokens=stop)
    for k, tid in enumerate(ids):
        if timings is not None:
            t0 = time.perf_counter_ns()
            m.fill_next_token_bitmask(bitmask)
            timings.append(time.perf_counter_ns() - t0)
        if not m.accept_token(tid):
            return False, k
    return m.is_terminated(), len(ids)


def llg_run(tok, grammar, ids, bitmask, timings=None):
    m = llg.LLMatcher(tok, grammar, log_level=0)
    for k, tid in enumerate(ids):
        if timings is not None:
            t0 = time.perf_counter_ns()
            lln.fill_next_token_bitmask(m, bitmask)
            timings.append(time.perf_counter_ns() - t0)
        if not m.consume_token(tid):
            return False, k
    return (m.is_stopped() or m.is_accepting()) and not m.is_error(), len(ids)


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


def main():
    path = sys.argv[1]
    n = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 10**9
    names = sys.argv[sys.argv.index("--models") + 1].split(",") if "--models" in sys.argv else list(MODELS)
    D = json.load(open(path))
    rng = random.Random(0)
    report = {}
    for name in names:
        repo, stop = MODELS[name]
        tk = AutoTokenizer.from_pretrained(repo)
        stop = stop or [tk.eos_token_id]
        V = max(len(tk), max(tk.get_vocab().values()) + 1)
        info = xgr.TokenizerInfo.from_huggingface(tk, vocab_size=V, stop_token_ids=stop)
        comp = xgr.GrammarCompiler(info, max_threads=1, cache_enabled=False)
        lt = llguidance.hf.from_tokenizer(tk, n_vocab=V, eos_token=stop[0])
        r = {"vocab": V, "stop": stop}
        for ds, pat in D["patterns"].items():
            t0 = time.perf_counter()
            cg = comp.compile_regex(pat)
            r[f"{ds}_xgr_compile_ms"] = round(1e3 * (time.perf_counter() - t0), 1)
            t0 = time.perf_counter()
            lgram = llg.LLMatcher.grammar_from_regex(pat)
            msg = llg.LLMatcher.validate_grammar(lgram, lt)
            ok_g = not msg
            llg.LLMatcher(lt, lgram, log_level=0)  # construction cost is part of the timing
            r[f"{ds}_llg_compile_ms"] = round(1e3 * (time.perf_counter() - t0), 1)
            r[f"{ds}_llg_valid"] = [ok_g, msg[:200]]
            texts = D["texts"][ds]
            if len(texts) > n:
                texts = rng.sample(texts, n)
            bm_x = xgr.allocate_token_bitmask(1, V)
            bm_l = lln.allocate_token_bitmask(1, V)
            acc_x = acc_l = 0
            fails = []
            ntok = []
            tx, tl = [], []
            for idx, t in enumerate(texts):
                ids = tk.encode(t, add_special_tokens=False) + [stop[0]]
                ntok.append(len(ids))
                time_it = idx < 200
                ox, _ = xgr_run(cg, info, ids, stop, bm_x, tx if time_it else None)
                ol, _ = llg_run(lt, lgram, ids, bm_l, tl if time_it else None)
                acc_x += ox
                acc_l += ol
                if not (ox and ol) and len(fails) < 5:
                    fails.append((ox, ol, t[:120]))
            r[f"{ds}_n"] = len(texts)
            r[f"{ds}_tokens_median"] = st.median(ntok)
            r[f"{ds}_accept_xgr"] = acc_x
            r[f"{ds}_accept_llg"] = acc_l
            r[f"{ds}_fails"] = fails
            r[f"{ds}_mask_us_xgr"] = {"mean": round(st.mean(tx) / 1e3, 1), "p50": round(q(tx, .5) / 1e3, 1),
                                      "p99": round(q(tx, .99) / 1e3, 1), "max": round(max(tx) / 1e3, 1), "n": len(tx)}
            r[f"{ds}_mask_us_llg"] = {"mean": round(st.mean(tl) / 1e3, 1), "p50": round(q(tl, .5) / 1e3, 1),
                                      "p99": round(q(tl, .99) / 1e3, 1), "max": round(max(tl) / 1e3, 1), "n": len(tl)}
            # malformed: must be rejected by both (not accepted as complete)
            bad = {}
            for t in texts[:50]:
                for kind, b in malformed(t).items():
                    if b == t:
                        continue
                    ids = tk.encode(b, add_special_tokens=False) + [stop[0]]
                    ox, _ = xgr_run(cg, info, ids, stop, bm_x)
                    ol, _ = llg_run(lt, lgram, ids, bm_l)
                    bad.setdefault(kind, [0, 0, 0])
                    bad[kind][0] += 1
                    bad[kind][1] += ox
                    bad[kind][2] += ol
            r[f"{ds}_malformed_accepted(n,xgr,llg)"] = bad
            # allowed-token counts along a path: how constraining is the grammar per step
            t = texts[0]
            ids = tk.encode(t, add_special_tokens=False) + [stop[0]]
            m = xgr.GrammarMatcher(cg, override_stop_tokens=stop)
            allowed = []
            for tid in ids:
                m.fill_next_token_bitmask(bm_x)
                bits = np.unpackbits(bm_x[0].numpy().view(np.uint8), bitorder="little")[:V]
                allowed.append(int(bits.sum()))
                m.accept_token(tid)
            r[f"{ds}_allowed_per_step"] = {"min": min(allowed), "median": st.median(allowed), "max": max(allowed),
                                           "first_12": allowed[:12]}
        report[name] = r
        print(json.dumps({name: r}, indent=1), flush=True)
    json.dump(report, open(path.replace("targets_texts.json", "grammar_report.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
