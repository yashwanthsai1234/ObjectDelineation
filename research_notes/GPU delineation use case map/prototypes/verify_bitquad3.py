import random
random.seed(1)
N4=[(1,0),(-1,0),(0,1),(0,-1)]; N8=N4+[(1,1),(1,-1),(-1,1),(-1,-1)]
def ncomp(g,val,nb):
    H,W=len(g),len(g[0]); seen=set(); n=0
    for i in range(H):
        for j in range(W):
            if g[i][j]==val and (i,j) not in seen:
                n+=1; st=[(i,j)]; seen.add((i,j))
                while st:
                    y,x=st.pop()
                    for dy,dx in nb:
                        yy,xx=y+dy,x+dx
                        if 0<=yy<H and 0<=xx<W and g[yy][xx]==val and (yy,xx) not in seen:
                            seen.add((yy,xx)); st.append((yy,xx))
    return n
def pad(b):
    W=len(b[0]); return [[0]*(W+2)]+[[0]+r+[0] for r in b]+[[0]*(W+2)]
ok4=ok8=okA=okP=okV=0;T=2000
for t in range(T):
    h,w=random.randint(1,10),random.randint(1,10); pr=random.uniform(.2,.8)
    b=[[1 if random.random()<pr else 0 for _ in range(w)] for _ in range(h)]
    p=pad(b); H,W=len(p),len(p[0])
    n=[0]*5; nd=0
    for i in range(H-1):
        for j in range(W-1):
            a,bb,c,d=p[i][j],p[i][j+1],p[i+1][j],p[i+1][j+1]; s=a+bb+c+d
            if s==2 and a==d: nd+=1
            else: n[s]+=1
    n1,n2,n3,n4=n[1],n[2],n[3],n[4]
    area=sum(map(sum,b))
    fg4,fg8=ncomp(p,1,N4),ncomp(p,1,N8)
    bg4,bg8=ncomp(p,0,N4)-1,ncomp(p,0,N8)-1   # holes = bounded bg comps
    if 4*(fg4-bg8)==n1-n3+2*nd: ok4+=1
    if 4*(fg8-bg4)==n1-n3-2*nd: ok8+=1
    if 4*area==n1+2*n2+3*n3+4*n4+2*nd: okA+=1
    edges=sum(abs(p[i][j]-p[i][j+1]) for i in range(H) for j in range(W-1))+sum(abs(p[i][j]-p[i+1][j]) for i in range(H-1) for j in range(W))
    if edges==n1+n2+n3+2*nd: okP+=1
    # rectilinear corners: count boundary grid vertices where boundary turns; each Q1/Q3 = 1 turn, diag = 2 turns
    # independent check: corners = vertices where the incident unit boundary edges are not collinear
    corners=0
    for i in range(H-1):
        for j in range(W-1):
            a,bb,c,d=p[i][j],p[i][j+1],p[i+1][j],p[i+1][j+1]
            top=a!=bb; bot=c!=d; left=a!=c; right=bb!=d   # 4 unit edges meeting at the quad centre vertex
            k=top+bot+left+right
            if k==2 and not ((top and bot) or (left and right)): corners+=1
            elif k==4: corners+=2
    if corners==n1+n3+2*nd: okV+=1
print(f"chi_fg4 (holes=8-conn bg) == (n1-n3+2nD)/4 : {ok4}/{T}")
print(f"chi_fg8 (holes=4-conn bg) == (n1-n3-2nD)/4 : {ok8}/{T}")
print(f"area == (n1+2n2+3n3+4n4+2nD)/4           : {okA}/{T}")
print(f"edge-count == n1+n2+n3+2nD (Gray P)      : {okP}/{T}")
print(f"turn-vertices == n1+n3+2nD                : {okV}/{T}")
