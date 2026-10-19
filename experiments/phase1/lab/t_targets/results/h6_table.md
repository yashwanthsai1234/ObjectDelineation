| masks | pixel-frame simplifier | S divides 1000 | S > 1000 | n | masks with flipped px % | flipped px / mask | IoU(rounded vs unrounded) min | topology changed % | became illegal % | identity round trip exact % |
|---|---|---|---|---|---|---|---|---|---|---|
| blob | dp1px | no | no | 84 | 100.0 | 133.8 | 0.9884 | 0.0 | 0.0 | - |
| blob | dp1px | no | yes | 24 | 100.0 | 1043.9 | 0.9931 | 0.0 | 0.0 | - |
| blob | dp1px | yes | no | 72 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | - |
| blob | identity | no | no | 84 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| blob | identity | no | yes | 24 | 100.0 | 516.0 | 0.9958 | 0.0 | 12.5 | 0.0 |
| blob | identity | yes | no | 72 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| blob | minseg1px | no | no | 84 | 100.0 | 123.7 | 0.9903 | 0.0 | 0.0 | - |
| blob | minseg1px | no | yes | 24 | 100.0 | 1029.6 | 0.9931 | 0.0 | 0.0 | - |
| blob | minseg1px | yes | no | 72 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | - |
| lidc | identity | no | no | 190 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| lidc | identity | yes | no | 113 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| lidc | minseg1px | no | no | 190 | 71.1 | 1.8 | 0.8738 | 0.0 | 0.0 | - |
| lidc | minseg1px | yes | no | 113 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | - |
| riga | dp1px | no | no | 78 | 100.0 | 46.9 | 0.9919 | 0.0 | 0.0 | - |
| riga | dp1px | no | yes | 40 | 100.0 | 120.1 | 0.9800 | 2.5 | 0.0 | - |
| riga | dp1px | yes | no | 68 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | - |
| riga | identity | no | no | 78 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| riga | identity | no | yes | 40 | 100.0 | 120.7 | 0.9845 | 0.0 | 20.0 | 0.0 |
| riga | identity | yes | no | 68 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | 100.0 |
| riga | minseg1px | no | no | 78 | 100.0 | 45.4 | 0.9938 | 0.0 | 0.0 | - |
| riga | minseg1px | no | yes | 40 | 100.0 | 121.0 | 0.9812 | 0.0 | 0.0 | - |
| riga | minseg1px | yes | no | 68 | 0.0 | 0.0 | 1.0000 | 0.0 | 0.0 | - |
