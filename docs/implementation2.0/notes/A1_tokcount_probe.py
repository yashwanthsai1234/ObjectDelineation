import json, glob, sys, numpy as np
from transformers import AutoTokenizer
tok = AutoTokenizer.from_pretrained(glob.glob(sys.argv[1])[0])
for ds in ("lidc","riga"):
    for sp in ("train","val"):
        try: recs=[json.loads(l) for l in open(f"data/prepared/{ds}/{sp}.jsonl")]
        except FileNotFoundError: continue
        L=[len(tok(t,add_special_tokens=False)["input_ids"])+1 for r in recs for t in r["targets"] if t]
        L=np.array(L)
        print(ds,sp,len(recs),"targets",len(L),"mean",L.mean().round(1),"p50",np.percentile(L,50),"p99",np.percentile(L,99),"p99.9",np.percentile(L,99.9),"max",L.max(),
              ">=600",(L>=600).sum(),">=800",(L>=800).sum(),">1000",(L>1000).sum(),">1200",(L>1200).sum(),">1536",(L>1536).sum())
