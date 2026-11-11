"""tools.profiling.c0_summary (plan C0 steps 5 and 6): U-MS4 (the decode-tail formula) and U-MS5 (the
c0_summary.json schema), on nsys_baseline's hand-made export (tests/goldens/nsys_small)."""

import json

import pytest

from tests.experiments.test_nsys_baseline import export, log, ranges
from tools.profiling import c0_summary as C0

LENGTHS = {"16": [1536, 1200, 130, 120, 110, 100, 90, 80, 70], "17": [100] * 8}


# U-MS4 -----------------------------------------------------------------------------------------------------------


def test_decode_tail_is_4_3_ms_per_token_of_longest_minus_8th_longest():
    assert C0.decode_tail_s([1536, 1200, 130, 120, 110, 100, 90, 80, 70]) == pytest.approx(0.0043 * (1536 - 80))
    assert C0.decode_tail_s([80, 1536, 90, 70, 1200, 100, 130, 110, 120]) == pytest.approx(6.2608)  # order-free
    assert C0.decode_tail_s([5] * 8) == 0.0
    assert C0.decode_tail_s(list(range(8))) == pytest.approx(0.0043 * 7)  # exactly 8: the 8th is the shortest


def test_decode_tail_needs_eight_answers():
    with pytest.raises(ValueError, match="7 answers"):
        C0.decode_tail_s([1536, 1200, 130, 120, 110, 100, 90])


# U-MS5 -----------------------------------------------------------------------------------------------------------


def run(tmp_path, csv=None, window=(16, 17), lengths=LENGTHS, logs=None):
    stats, db = export(tmp_path, csv)
    return C0.summary(stats, db, logs or log([10.0] * 17), window, lengths, r"^k2$")


def test_hand_values(tmp_path):
    logs = log([10.0] * 17)
    logs[15] |= {"perf/host_syncs/batch": 3, "perf/host_syncs/step": 9, "perf/host_syncs_vllm": 2}
    doc = run(tmp_path, logs=logs)
    s1, s2 = doc["steps"]
    assert doc["window"] == [16, 17] and doc["graph_tables"] == ["CUPTI_ACTIVITY_KIND_GRAPH_TRACE"]
    assert s1["wall_s"] == pytest.approx({"step": 100e-9, "generate": 40e-9, "batch": 25e-9, "sync": 15e-9,
                                          "rewards": 30e-9})  # fmt: skip
    assert s1["gpu_active"] == pytest.approx({"step": 0.55, "generate": 0.75, "batch": 0.8, "sync": 10 / 15,
                                              "rewards": 1 / 3})  # fmt: skip
    assert s2["gpu_active"] == pytest.approx({"step": 0.4, "forward": 0.4, "generate": 0.5})  # the graph run counts
    assert s1["syncs_by_range"] == {"step": 6, "generate": 1, "batch": 1, "sync": 1, "rewards": 1}
    assert s2["syncs_by_range"] == {"step": 1, "forward": 1, "generate": 1}  # each range on its own thread
    assert s1["syncs"] == {"torch": {"batch": 3, "step": 9}, "torch_vllm": 2, "nsys_trainer_thread": 4, "nsys_vllm": 2}
    assert s2["syncs"] == {"torch": {}, "torch_vllm": None, "nsys_trainer_thread": 1, "nsys_vllm": 1}
    assert s1["h2d_bytes"] == {"trainer": 1000, "vllm": 0} and s1["d2h_bytes"] == {"trainer": 0, "vllm": 0}
    assert s1["input_processing_s"] == pytest.approx(20e-9) and s2["input_processing_s"] is None  # 6a
    assert s1["decode_tail_s"] == pytest.approx(6.2608) and s2["decode_tail_s"] == 0.0  # 6b
    assert s1["generate_minus_handovers_s"] == pytest.approx(40e-9) and s1["longest_answer"] == 1536  # 6c
    assert s2["longest_answer"] == 100


def test_handovers_reward_thread_and_copy_directions(tmp_path):
    """A sync and a wake inside generate leave it; a sync on another thread (the reward thread) is in neither nsys
    count of the trainer thread; DtoH inside generate is vLLM's, HtoD outside is the trainer's."""
    gpu = ["40,5,1,,,,,,,,,,0.002,1.0,Device,Pageable,dev,1,,7,[CUDA memcpy DtoH]",
           "70,5,1,,,,,,,,,,0.5,1.0,Pageable,Device,dev,1,,7,[CUDA memcpy Host-to-Device]"]  # fmt: skip
    csv = {"nvtx_pushpop_trace": ranges((0, 100, "step", 1), (0, 50, "generate", 1), (10, 20, "sync", 1),
                                        (30, 35, "wake", 1)),  # fmt: skip
           "cuda_api_trace": ["60,1,cudaStreamSynchronize,0,1,7,3,20,reward",
                              "15,1,cudaDeviceSynchronize,0,2,7,1,20,python"],  # fmt: skip
           "nvtx_gpu_proj_trace": [], "cuda_gpu_trace": gpu}
    (s,) = run(tmp_path, csv=csv, window=(16, 16), lengths={"16": list(range(8))})["steps"]
    assert s["generate_minus_handovers_s"] == pytest.approx(35e-9)
    assert s["syncs"]["nsys_trainer_thread"] == 0 and s["syncs"]["nsys_vllm"] == 1
    assert s["h2d_bytes"] == {"trainer": 500_000, "vllm": 0} and s["d2h_bytes"] == {"trainer": 0, "vllm": 2000}
    assert s["decode_tail_s"] == pytest.approx(0.0043 * 7)


@pytest.mark.parametrize("field", [*C0.FIELDS, *(f"syncs.{f}" for f in C0.SYNC_FIELDS)])
def test_a_missing_field_fails(tmp_path, field):
    doc = run(tmp_path)
    C0.check(doc)
    head, _, sub = field.partition(".")
    row = doc["steps"][1]
    del (row[head] if sub else row)[sub or head]
    with pytest.raises(ValueError, match=rf"steps\[1\]: {field}"):
        C0.check(doc)


@pytest.mark.parametrize("key", ["window", "graph_tables"])
def test_a_missing_top_level_field_fails(tmp_path, key):
    doc = run(tmp_path)
    del doc[key]
    with pytest.raises(ValueError, match=key):
        C0.check(doc)


def test_main_writes_the_summary_and_refuses_missing_lengths(tmp_path):
    stats, db = export(tmp_path)
    (tmp_path / "rl.json").write_text(json.dumps({"log": log([10.0] * 17)}))
    (tmp_path / "len.json").write_text(json.dumps(LENGTHS))
    argv = ["--stats", stats, "--sqlite", str(db), "--log", str(tmp_path / "rl.json"), "--window", "16", "17",
            "--lengths", str(tmp_path / "len.json"), "--vision-kernel", "k2", "--out", str(tmp_path / "c0.json")]
    assert C0.main(argv) == 0
    doc = json.loads((tmp_path / "c0.json").read_text())
    assert [s["step"] for s in doc["steps"]] == [16, 17] and doc["steps"][0]["decode_tail_s"] == pytest.approx(6.2608)
    (tmp_path / "len.json").write_text(json.dumps({"16": LENGTHS["16"]}))
    with pytest.raises(SystemExit, match="no answer lengths for step 17"):
        C0.main(argv)
    (tmp_path / "len.json").write_text(json.dumps({"16": LENGTHS["16"], "17": [1] * 7}))
    with pytest.raises(SystemExit, match="7 answers"):
        C0.main(argv)
