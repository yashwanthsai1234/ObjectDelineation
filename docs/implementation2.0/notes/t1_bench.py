import glob,time,hashlib,os,io
import numpy as np, pyarrow as pa, pyarrow.parquet as pq
from PIL import Image
fs=sorted(glob.glob('data/views/lidc-train-shifted/images/lidc/*/*/*.png'))[:300]
g=[np.asarray(Image.open(f))[...,0] for f in fs]
raw=g[0].nbytes
for lvl in (1,3,6):
    cs=[]
    for x in g[:100]:
        b=io.BytesIO();pq.write_table(pa.table({'a':pa.array([x.tobytes()],pa.binary())}),b,compression='zstd',compression_level=lvl,use_dictionary=False);cs.append(b.getvalue())
    t=time.time();[pq.read_table(io.BytesIO(c)) for c in cs];dt=(time.time()-t)/100
    print('zstd',lvl,'gray bytes',np.mean([len(c) for c in cs]),'dec ms',dt*1e3)
t=time.time();[hashlib.blake2b(x.tobytes()).hexdigest() for x in g];print('blake2b gray ms',(time.time()-t)/300*1e3)
t=time.time();[hashlib.sha256(x.tobytes()).hexdigest() for x in g];print('sha256 gray ms',(time.time()-t)/300*1e3)
t=time.time();[np.asarray(Image.open(f)) for f in fs];print('png ms',(time.time()-t)/300*1e3)
print('png mean B',np.mean([os.path.getsize(f) for f in fs]),'gray raw',raw)
print('22023 views gray GB',22023*raw/1e9,'RGB GB',22023*raw*3/1e9, 'rows 24624 gray GB',24624*raw/1e9)
print('png total shifted+rescaled GB',14682*162771/1e9)
