from __future__ import print_function, division
from collections import Counter
from frisk.kmerhash import *
import screed
from matplotlib import pyplot as plt
from functools import partial
import multiprocessing as mp
from os import path


def seq_ivom(gen):
    kv = build_kmer_vec(1, K)
    hash_seq(str(gen), kv)
    return ivom(kv)


def window_kli(start, seq, giv, winsz):
    winseq = seq[start:start+winsz]
    iv = seq_ivom(winseq)
    return start, kli(giv, iv)

if __name__ == "__main__":
    import sys
    K=8
    winsz = 5000
    n = 10000000

    r = screed.open(sys.argv[1])
    for x in r:
        seq = x.sequence[:n]
        break

    if path.exists('genome.npy'):
        giv = np.load('genome.npy')
        print("load genome IVOM")
    else:
        giv = seq_ivom(seq)
        np.save('genome.npy', giv)
        print("calculated genome IVOM")

    do_window = partial(window_kli, seq=seq, giv=giv, winsz=winsz)

    klis = []
    print(len(seq))
    starts = range(0, len(seq) - winsz + 1, int(winsz/2))
    if True:
        pool = mp.Pool()
        for window in pool.imap(do_window, starts):
            klis.append(window)
            if len(klis) % 100 == 0:
                print("done", len(klis), "windows")
        klis = np.array([k for _, k in sorted(klis)])
    else:
        for start in starts:
            print(start)
            winseq = seq[start:start+winsz]
            c = Counter()
            for base in winseq:
                c[base]+=1
            for  k, v in c.most_common():
                print(k, v/len(winseq))
            iv = seq_ivom(winseq)
            print(sum(iv), len(winseq))
            klis.append(kli(giv, iv))
            print(klis[-1])
            if len(klis) % 100 == 0:
                print("done", len(klis), "windows")


    #np.save('window_klis.npy', klis)

    plt.plot(klis)
    plt.savefig('frisk.pdf')
