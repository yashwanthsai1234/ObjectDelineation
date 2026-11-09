import re, sys
T = open(sys.argv[1]).read()
C = dict(svc=("rgba(6,78,59,0.4)","#34d399"), db=("rgba(76,29,149,0.4)","#a78bfa"), q=("rgba(251,146,60,0.3)","#fb923c"),
         ext=("rgba(30,41,59,0.5)","#94a3b8"), sec=("rgba(136,19,55,0.4)","#fb7185"), cache=("rgba(8,51,68,0.4)","#38bdf8"),
         cloud=("rgba(120,53,15,0.3)","#fbbf24"))
def box(x,y,w,h,kind,title,*lines):
    f,s = C[kind]; o=[f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{f}" stroke="{s}" stroke-width="1.5"/>',
        f'<text x="{x+w/2}" y="{y+17}" fill="white" font-size="11" font-weight="600" text-anchor="middle">{title}</text>']
    for i,l in enumerate(lines): o.append(f'<text x="{x+w/2}" y="{y+32+i*12}" fill="#94a3b8" font-size="8.5" text-anchor="middle">{l}</text>')
    return "\n".join(o)
def arr(x1,y1,x2,y2,col="#64748b",label=None,dash=False,lx=None,ly=None):
    d=' stroke-dasharray="5,4"' if dash else ''
    o=f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" stroke-width="1.4"{d} marker-end="url(#arrowhead)"/>'
    if label: o+=f'\n<text x="{lx if lx is not None else (x1+x2)/2+4}" y="{ly if ly is not None else (y1+y2)/2-4}" fill="#94a3b8" font-size="8">{label}</text>'
    return o
P=[]
# boundaries
P.append('<rect x="20" y="40" width="560" height="620" rx="12" fill="rgba(52,211,153,0.03)" stroke="#34d399" stroke-width="1" stroke-dasharray="8,4"/>')
P.append('<text x="34" y="58" fill="#34d399" font-size="10" font-weight="600">MAC  (CPU, $0, gpucheck CPU build)</text>')
P.append('<rect x="610" y="40" width="370" height="620" rx="12" fill="rgba(251,191,36,0.04)" stroke="#fbbf24" stroke-width="1" stroke-dasharray="8,4"/>')
P.append('<text x="624" y="58" fill="#fbbf24" font-size="10" font-weight="600">VAST BOX  (RTX 5090, nvcc build)</text>')
# arrows first
P+= [arr(150,120,150,148), arr(150,200,150,228,label="certify",lx=156),
     arr(210,250,330,250,label="draw",lx=255), arr(390,272,390,298,label="SFT",lx=396),
     arr(150,272,150,358), arr(390,350,390,378), arr(330,404,250,404,lx=270),
     arr(250,175,330,175,dash=True,label="flat blob",lx=262),
     arr(150,410,150,458), arr(390,430,300,470,lx=330),
     arr(220,488,330,560,label="pass",lx=262,ly=520), arr(150,516,150,560,col="#fb7185",dash=True,label="fail = no version",lx=156),
     arr(450,585,628,585,col="#fbbf24",label="hf download (pinned rev)",lx=470,ly=578),
     arr(720,608,720,525+0,col="#fbbf24"),
     arr(720,470,720,440), arr(690,380,690,350,label="32 spec ids",lx=696),
     arr(780,320,840,320,dash=True), arr(700,290,700,250), arr(860,290,860,250),
     arr(700,200,700,170), arr(760,140,820,140), arr(860,200,860,170),
     arr(450,180,628,320,col="#38bdf8",dash=True,label="contexts 476 MB",lx=500,ly=232)]
# boxes
P+= [box(80,80,140,40,"ext","RAW (untrusted)","RIGA 26 GB, LIDC 5.2 GB"),
     box(80,148,170,52,"db","L1 ITEM STORE","Parquet, 5.08 GB, SHA-256","+ nodule_id"),
     box(80,228,130,44,"db","L2 DOCTOR RECORDS","bits 2 MB + invariants"),
     box(330,148,120,52,"db","CONTEXT BLOB","int16 HU, 476 MB","+ index"),
     box(330,228,120,44,"db","L3 VIEW SPECS","~2 MB, spec_id key"),
     box(330,298,120,52,"svc","TARGETS (SFT)","CUDA rings -> simplify","-> certify"),
     box(330,378,120,52,"db","FROZEN EXAMS","PNG + caches","exam_version"),
     box(80,358,140,52,"db","L4 ROW SETS","RL: ids + weights","SFT: spec + target"),
     box(80,458,140,58,"sec","AUDIT G1-G15","coverage vs VAL","leakage, parity refs"),
     box(250,560,200,50,"q","data.lock.json + datacard","version = sha256(parents,","code, gate results)"),
     box(640,470,160,52,"sec","PREFLIGHT G1 G8 G9","hashes, CPU==nvcc tables","render bytes"),
     box(640,608,120,40,"cloud","HF private repo","LIDC only, Xet dedup"),
     box(640,380,170,52,"svc","SAMPLER","cell-matched weights","logs spec ids"),
     box(640,290,140,60,"cache","PINNED HOST CTX","476 MB, H2D ~1 MB/step"),
     box(800,290,160,60,"cache","DEVICE DOCTORS","~2 MB resident","tables built per step"),
     box(640,200,120,50,"svc","RENDER kernel","window, op, bicubic","uint8 512^2"),
     box(800,200,160,50,"svc","TABLES kernel","place bits, op, EDT","freed after step"),
     box(640,110,120,60,"svc","GPU pixel build","(Phase C2) -> vLLM","answers"),
     box(820,110,140,60,"svc","SCORER (gpucheck)","IoU vs doctor set","rewards + data_version")]
# legend
lg=[("svc","service / kernel"),("db","store"),("cache","resident / pinned"),("sec","gate"),("q","version lock"),("cloud","transport")]
for i,(k,t) in enumerate(lg):
    f,s=C[k]; y=560+i*16
    P.append(f'<rect x="830" y="{y}" width="14" height="10" rx="2" fill="{f}" stroke="{s}" stroke-width="1"/><text x="850" y="{y+9}" fill="#94a3b8" font-size="8.5">{t}</text>')
svg_body="\n".join(P)
T = re.sub(r'(<rect width="100%" height="100%" fill="url\(#grid\)" />).*?(</svg>)', lambda m: m.group(1)+"\n"+svg_body+"\n"+m.group(2), T, count=1, flags=re.S)
T = T.replace("[PROJECT NAME] Architecture","objdel data system (SD3, data-centric)").replace("[Subtitle description]",
    "Raw LIDC/RIGA -> certified doctor records -> view specs -> gated, versioned row sets -> GPU-fed training store (up to Phase E)")
T = T.replace('<svg viewBox="0 0 1000 680">','<svg viewBox="0 0 1000 680">',1)
cards=[("rose","Gates that block a version",["G10 coverage vs val: min cell ratio >= 0.8 (today 0.48)","G3/G4 leakage = 0 (patient, series, pixels)","G8/G9 CPU==CUDA tables, targets, render bytes","G15 every result names data_version"]),
       ("amber","Key decisions",["Specs (KB) are the dataset; pixels derived on box","Doctor tables built per step from 2 MB of bits","Contexts pinned on host, not in VRAM","Exams frozen with their own exam_version"]),
       ("violet","Weakest point / next",["A2 hypothesis: window tables = f(record, offset)","Proven only by G8 on every window","3D phase: contexts ~3.7 GB, streaming needed","No GPU decode until compressed DICOM volumes"])]
html=""
for dot,title,items in cards:
    html+=f'<div class="card"><div class="card-header"><div class="card-dot {dot}"></div><h3>{title}</h3></div><ul>'+"".join(f"<li>&bull; {i}</li>" for i in items)+"</ul></div>\n"
T = re.sub(r'(<div class="cards">).*?(</div>\s*<!-- Footer -->)', lambda m: m.group(1)+"\n"+html+"</div>\n\n    <!-- Footer -->", T, count=1, flags=re.S)
T = re.sub(r'\[Project Name\] &bull; \[Additional metadata\]|\[Project Name\] • \[Additional metadata\]', "SD3_data.md &bull; 2026-11-08 &bull; numbers tagged in the doc", T)
open(sys.argv[2],"w").write(T)
