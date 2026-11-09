"""Replace the text between two exact marker lines (start kept, end kept) of a file with the text of another file."""
import sys
path, start, end, newf = sys.argv[1:5]
L = open(path).read().split("\n")
i = L.index(start); j = L.index(end, i + 1)
new = open(newf).read().rstrip("\n").split("\n")
open(path, "w").write("\n".join(L[: i + 1] + new + L[j:]))
print("replaced", j - i - 1, "lines with", len(new))
