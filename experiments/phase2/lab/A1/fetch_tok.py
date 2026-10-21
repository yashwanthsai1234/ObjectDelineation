import sys

from huggingface_hub import hf_hub_download, list_repo_files

repos = sys.argv[1:]
for r in repos:
    try:
        files = list_repo_files(r)
    except Exception as e:
        print(r, "LIST-FAIL", type(e).__name__, str(e)[:100])
        continue
    want = [
        f
        for f in files
        if f.split("/")[-1]
        in (
            "tokenizer.json",
            "tokenizer_config.json",
            "vocab.json",
            "merges.txt",
            "special_tokens_map.json",
            "added_tokens.json",
            "tokenizer.model",
            "config.json",
            "preprocessor_config.json",
            "processor_config.json",
            "chat_template.jinja",
            "chat_template.json",
        )
    ]
    ok = []
    for f in want:
        try:
            hf_hub_download(r, f)
            ok.append(f)
        except Exception as e:
            print(r, f, "FAIL", type(e).__name__)
    print(r, "OK", ok)
