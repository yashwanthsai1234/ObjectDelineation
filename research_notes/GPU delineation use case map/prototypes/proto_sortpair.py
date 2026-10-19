# CPU reference of the proposed data-parallel design: every step is a map, sort, scan, or pointer-jump.
import random
random.seed(7)
N4=[(1,0),(-1,0),(0,1),(0,-1)]; N8=N4+[(1,1),(1,-1),(-1,1),(-1,-1)]
def label(g,val,nb):
    H,W=len(g),len(g[0]); L=[[0]*W for _ in range(H)]; n=0
    for i in range(H):
        for j in range(W):
            if g[i][j]==val and not L[i][j]:
                n+=1; st=[(i,j)]; L[i][j]=n
                while st:
                    y,x=st.pop()
                    for dy,dx in nb:
                        yy,xx=y+dy,x+dx
                        if 0<=yy<H and 0<=xx<W and g[yy][xx]==val and not L[yy][xx]:
                            L[yy][xx]=n; st.append((yy,xx))
    return L,n
def run(b,fg8):
    W0=len(b[0]); p=[[0]*(W0+2)]+[[0]+r+[0] for r in b]+[[0]*(W0+2)]; H,W=len(p),len(p[0])
    # 1) map: emit corners (x,y,h,v); h in {-1 L,+1 R}, v in {-1 U,+1 D}
    C=[]
    for i in range(H-1):
        for j in range(W-1):
            a,bq,c,d=p[i][j],p[i][j+1],p[i+1][j],p[i+1][j+1]
            up,dn,lf,rt=a!=bq,c!=d,a!=c,bq!=d
            x,y=j+1,i+1
            if up and dn and lf and rt:      # diagonal quad -> two corners
                diagA = (a==1)               # tl & br foreground
                # fg8: corners wrap the background pixels; fg4: wrap the foreground pixels
                wrap_tl_br = (diagA != fg8)  # which diagonal pair of pixels each corner wraps
                if wrap_tl_br: C += [(x,y,-1,-1),(x,y,+1,+1)]   # (L,U) around tl, (R,D) around br
                else:          C += [(x,y,+1,-1),(x,y,-1,+1)]   # (R,U) around tr, (L,D) around bl
            elif (up or dn) and (lf or rt) and (up+dn+lf+rt)==2:
                C.append((x,y,-1 if lf else +1,-1 if up else +1))
    n=len(C)
    # 2) sort by (y, x, h) and pair consecutive -> horizontal partner
    hp=[0]*n; o=sorted(range(n),key=lambda k:(C[k][1],C[k][0],C[k][2]))
    for t in range(0,n,2): u,v=o[t],o[t+1]; assert C[u][1]==C[v][1] and C[u][2]==+1 and C[v][2]==-1; hp[u],hp[v]=v,u
    # 3) sort by (x, y, v) and pair -> vertical partner
    vp=[0]*n; o=sorted(range(n),key=lambda k:(C[k][0],C[k][1],C[k][3]))
    for t in range(0,n,2): u,v=o[t],o[t+1]; assert C[u][0]==C[v][0] and C[u][3]==+1 and C[v][3]==-1; vp[u],vp[v]=v,u
    # 4) ring id via pointer jumping on the 2-regular graph: succ = vp[hp[k]] visits every other corner; min-label propagation
    succ=[vp[hp[k]] for k in range(n)]
    rid=list(range(n)); s=succ[:]
    for _ in range(64):                       # Wyllie-style: rid=min over orbit, log rounds
        rid=[min(rid[k],rid[s[k]],rid[hp[k]]) for k in range(n)]; s=[s[s[k]] for k in range(n)]
    for k in range(n): rid[k]=min(rid[k],rid[hp[k]])
    rings={}
    for k in range(n): rings.setdefault(rid[k],[]).append(k)
    # 5) classify by lexicographic-min corner: inside pixel (br of that vertex) fg => outer
    Lf,nf=label(p,1,N8 if fg8 else N4); Lb,nb=label(p,0,N4 if fg8 else N8)
    outer=holes=0; ok_parent=True
    for r,ks in rings.items():
        k=min(ks,key=lambda k:(C[k][1],C[k][0])); x,y=C[k][0],C[k][1]
        if p[y][x]==1: outer+=1
        else:
            holes+=1
            # parent = fg label of pixel above the hole's top-left-most corner (p[y-1][x])
            par=Lf[y-1][x]
            if par==0: par=Lf[y][x-1]
            # validate: with only the parent component as fg, the hole pixel is still enclosed
            g=[[1 if Lf[i][j]==par else 0 for j in range(W)] for i in range(H)]
            Lg,_=label(g,0,N4 if fg8 else N8)
            if Lg[y][x]==Lg[0][0]: ok_parent=False
    return outer==nf, holes==nb-1, ok_parent, n
T=1500; tot=[0,0,0]
for fg8 in (True,False):
    good=[0,0,0]
    for t in range(T):
        h,w=random.randint(1,14),random.randint(1,14); pr=random.uniform(.15,.85)
        b=[[1 if random.random()<pr else 0 for _ in range(w)] for _ in range(h)]
        r=run(b,fg8)
        for q in range(3): good[q]+=r[q]
    print(f"fg{'8' if fg8 else '4'}: outer==#fgCCL {good[0]}/{T}; holes==#bounded bgCCL {good[1]}/{T}; hole-parent valid {good[2]}/{T}")
