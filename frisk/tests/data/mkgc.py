from numpy import random


def randseq(length=1000000, gc=0.5):
    at = 1-gc
    p = [at/2, gc/2, gc/2, at/2]
    return "".join(list(random.choice(list("ACGT"), size=10000, p=p)))

seq = randseq(gc=0.2)
seq += randseq(length=10000, gc=0.8)
seq += randseq(gc=0.2)
with open("gc20gc80gc20.fa", "w") as fh:
    print(">gc20gc80gc20", seq, file=fh, sep='\n')
