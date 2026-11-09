import re,sys
for f,pat in [("B.md",r"\b(U\d+b?|P\d+|I\d+b?|E\d\d|G\d+[ab]?|F\d+|R\d+)\b"),("C.md",r"\b([UPIG]-[A-Z]{2,3}\d+[a-z]?|G-[A-Z]+\d*)\b"),("E.md",r"\b(U\d+|P\d+|I\d+)\b")]:
    t=open(f).read()
    # cut at review log
    body=t.split("Review log")[0] if False else t
    defined=set(re.findall(r"^\| ("+pat[2:-2].replace('(','(?:',1)+r") \|",body,re.M))
    defined=set(m for m in re.findall(r"^\|\s*\**([A-Z][A-Z0-9\-]*\d+[a-z]?)\**\s*\|",body,re.M))
    print(f,len(defined))
    refs=set(re.findall(pat,body))
    missing=sorted(r for r in refs if r not in defined)
    print(" missing:",missing[:200])
