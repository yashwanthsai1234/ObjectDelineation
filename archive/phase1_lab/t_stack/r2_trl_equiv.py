"""R2: the reward-v2 aggregator inside the exact TRL GRPOTrainer path equals rewardcore.batch.advantages_v2.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I r2_trl_equiv.py

B = 4 prompts x G = 8, hand-set term vectors injected through a provider (completions are from a tiny dummy LM
and ignored). TRL config: weight-0 per-term loggers + one aggregator (weight 1), multi_objective_aggregation=
"sum_then_normalize", scale_rewards="none", num_generations=8. The trainer's advantages are captured from
_generate_and_score_completions (and from trainer._logs["advantages"]) and compared with:
  run A (no drop): rewardcore.batch.advantages_v2 on a TERM_FIELDS matrix  -> |diff| <= 1e-6
  run B (drop_spread = 1e-2): low-spread + identical groups -> exactly 0, other groups == advantages_v2 on the
                              kept groups only.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("HF_HOME", os.path.join(os.path.dirname(os.path.abspath(__file__)), "hf"))
import numpy as np  # noqa: E402
import torch  # noqa: E402
from datasets import Dataset  # noqa: E402

import prod  # noqa: E402
from rewardcore.batch import TERM_FIELDS, advantages_v2  # noqa: E402

G = 8
SIG_A, SIG_B = prod.sig_id((0,)), prod.sig_id((1,))


def groups(variant):
    """Hand-set per-rollout terms for 4 groups: dict name -> (G,) arrays."""
    def g(iou, f5, legal, sig, freq, tany):
        return {"iou_tc": np.array(iou, float), "f5": np.array(f5, float), "legal": np.array(legal, float),
                "sig_id": np.array(sig, float), "topo_freq": np.array(freq, float), "topo_any": np.array(tany, float),
                "format_ok": np.ones(G)}
    illegal = g([.91, .88, .93, .90, .86, .92, .89, .95], [.97, .95, .98, .96, .93, .97, .94, .99],
                [1, 1, 1, 0, 1, 1, 1, 1], [SIG_A] * G, [1.0] * G, [1] * G)
    identical = g([.9] * G, [.95] * G, [1] * G, [SIG_A] * G, [1.0] * G, [1] * G)
    minority = g([.85, .87, .83, .88, .86, .84, .80, .89], [.93, .94, .92, .95, .93, .92, .90, .96],
                 [1] * G, [SIG_A] * 6 + [SIG_B] + [SIG_A], [2 / 3] * 6 + [1 / 3] + [2 / 3], [1] * G)
    if variant == "A":
        last = g([.70, .55, .82, .40, .65, .77, .30, .60], [.80, .70, .90, .50, .75, .85, .45, .72],
                 [1] * G, [SIG_A, SIG_A, SIG_A, SIG_B, SIG_A, SIG_A, SIG_B, SIG_A],
                 [1, 1, 1, 0, 1, 1, 0, 1], [1, 1, 1, 0, 1, 1, 0, 1])  # normal group (2 rollouts no rater drew)
    else:
        last = g([.900, .901, .903, .902, .904, .900, .905, .902], [.950, .951, .952, .950, .953, .951, .950, .952],
                 [1] * G, [SIG_A] * G, [1.0] * G, [1] * G)  # low spread (< 1e-2), one topology
    return [illegal, identical, minority, last]


def term_matrix(T):
    """Flat provider terms -> rewardcore TERM_FIELDS matrix (only the fields advantages_v2 reads matter)."""
    n = len(T["iou_tc"])
    M = np.zeros((n, len(TERM_FIELDS)))
    f = {k: i for i, k in enumerate(TERM_FIELDS)}
    M[:, f["iou_tc"]], M[:, f["band"]] = T["iou_tc"], T["f5"]
    M[:, f["topo_any"]], M[:, f["topo_freq"]], M[:, f["sig_id"]] = T["topo_any"], T["topo_freq"], T["sig_id"]
    return M


def run(variant, drop):
    from transformers import AutoModelForCausalLM, AutoTokenizer
    from trl import GRPOConfig, GRPOTrainer

    hand = groups(variant)
    seen = []

    def provider(prompts, completions, completion_ids, gid, **kw):
        gid = list(gid)
        blocks = [gid[i:i + G] for i in range(0, len(gid), G)]
        assert all(len(set(b)) == 1 for b in blocks), "groups are not contiguous"
        seen.append(gid)
        out = {k: np.concatenate([hand[b[0]][k] for b in blocks]) for k in hand[0]}
        return out

    agg = prod.TRLReward(provider, group=G, drop_spread=drop)
    funcs, weights = agg.funcs()
    mid = "trl-internal-testing/tiny-Qwen2ForCausalLM-2.5"
    tok = AutoTokenizer.from_pretrained(mid, padding_side="left")
    model = AutoModelForCausalLM.from_pretrained(mid, dtype=torch.float32)
    ds = Dataset.from_list([{"prompt": f"prompt {i}", "gid": i} for i in range(4)])
    cfg = GRPOConfig(output_dir=str(prod.HERE / "runs" / f"r2_{variant}"), max_steps=1, per_device_train_batch_size=32,
                     gradient_accumulation_steps=1, num_generations=G, max_completion_length=8, temperature=1.0,
                     multi_objective_aggregation="sum_then_normalize", scale_rewards="none", loss_type="dr_grpo",
                     reward_weights=weights, beta=0.0, use_cpu=True, bf16=False, report_to="none", logging_steps=1,
                     save_strategy="no", seed=0, gradient_checkpointing=False, remove_unused_columns=False)
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=cfg, train_dataset=ds, processing_class=tok)
    captured = []
    orig = trainer._generate_and_score_completions

    def hook(*a, **k):
        out = orig(*a, **k)
        captured.append(out["advantages"].detach().cpu().double().numpy().copy())
        return out
    trainer._generate_and_score_completions = hook
    logs = []
    from transformers import TrainerCallback

    class Sink(TrainerCallback):
        def on_log(self, args, state, control, logs_=None, **kw):
            pass
    trainer.train()
    logs = [h for h in trainer.state.log_history if "loss" in h]

    T = agg.history[-1]
    trl_adv = captured[-1]
    trl_logged_adv = np.array(list(trainer._logs["advantages"])[-len(trl_adv):], float)
    gid_order = [b[0] for b in [seen[-1][i:i + G] for i in range(0, len(seen[-1]), G)]]
    res = {"variant": variant, "drop_spread": drop, "group_order_seen": gid_order}
    if drop is None:
        ref = advantages_v2(term_matrix(T), group=G, c_topo=2.0, lam=0.5, legal=T["legal"])
    else:
        ref = prod.advantages(T["iou_tc"], T["f5"], T["legal"], T["sig_id"], T["topo_freq"], T["topo_any"],
                              group=G, drop_spread=drop)
        # kept groups must equal advantages_v2 run on the kept groups only
        keep = [i for i, b in enumerate(np.asarray(ref).reshape(-1, G)) if np.any(b != 0)]
        idx = np.concatenate([np.arange(i * G, (i + 1) * G) for i in keep])
        sub = {k: v[idx] for k, v in T.items() if isinstance(v, np.ndarray) and v.shape[0] == len(ref)}
        ref_kept = advantages_v2(term_matrix(sub), group=G, legal=sub["legal"])
        res["kept_groups_vs_advantages_v2_maxdiff"] = float(np.abs(np.asarray(ref)[idx] - ref_kept).max())
        nodrop = advantages_v2(term_matrix(T), group=G, legal=T["legal"]).reshape(-1, G)
        res["without_drop_lowspread_group_max_abs_adv"] = float(np.abs(nodrop[[gid_order.index(3)]]).max())
    ref = np.asarray(ref)
    res["max_abs_diff_trl_vs_ref"] = float(np.abs(trl_adv - ref).max())
    res["max_abs_diff_logged_vs_ref"] = float(np.abs(trl_logged_adv - ref).max())
    per_group = ref.reshape(-1, G)
    trl_g = trl_adv.reshape(-1, G)
    res["per_group"] = {f"g{gid_order[i]}": {"ref": per_group[i].round(4).tolist(), "trl": trl_g[i].round(4).tolist()}
                        for i in range(len(gid_order))}
    ident = gid_order.index(1)
    res["identical_group_trl_max_abs"] = float(np.abs(trl_g[ident]).max())
    if drop is not None:
        res["lowspread_group_trl_max_abs"] = float(np.abs(trl_g[gid_order.index(3)]).max())
    ill = gid_order.index(0)
    res["illegal_rollout_adv"] = float(trl_g[ill][3])
    res["illegal_is_min_in_group"] = bool(trl_g[ill][3] == trl_g[ill].min())
    mino = gid_order.index(2)
    res["minority_rollout_adv_vs_group_mean_of_others"] = float(trl_g[mino][6] - np.delete(trl_g[mino], 6).mean())
    # per-term logs == the hand-set vectors (in TRL order)
    term_ok, mean_ok = {}, {}
    last = logs[-1]
    for name in prod.LOG_TERMS:
        logged = np.array(list(trainer._logs["rewards"][name])[-len(ref):], float)
        term_ok[name] = float(np.abs(logged - T[name]).max())
        mean_ok[name] = float(abs(last[f"rewards/{name}/mean"] - np.mean(T[name])))
    res["per_term_log_max_abs_diff"] = term_ok
    res["per_term_mean_log_abs_diff"] = mean_ok
    res["aggregate_logged_mean"] = last.get("rewards/v2_aggregate/mean")
    res["logged_reward"] = last.get("reward")
    res["nan_in_logs"] = any(isinstance(v, float) and np.isnan(v) for v in last.values())
    res["PASS"] = bool(res["max_abs_diff_trl_vs_ref"] <= 1e-6 and max(term_ok.values()) <= 1e-6
                       and max(mean_ok.values()) <= 1e-6 and res["identical_group_trl_max_abs"] <= 1e-6
                       and (drop is None or (res["lowspread_group_trl_max_abs"] == 0.0
                                             and res["kept_groups_vs_advantages_v2_maxdiff"] <= 1e-9)))
    print(json.dumps(res, indent=1), flush=True)
    return res


if __name__ == "__main__":
    out = {"A_no_drop": run("A", None), "B_drop": run("B", 1e-2)}
    (prod.HERE / "results").mkdir(exist_ok=True)
    (prod.HERE / "results" / "r2_trl_equiv.json").write_text(json.dumps(out, indent=1))
    print("R2", "PASS" if all(v["PASS"] for v in out.values()) else "FAIL")
