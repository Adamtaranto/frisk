#!/usr/bin/env python
#Pybedtools.py

from operator import itemgetter
import pybedtools

infile = "Parastagonospora_nodorum_SN15.gff3"


def gffFilter(rec, **kwargs):
	if 'feature' in kwargs.keys():
		if rec[2] in kwargs['feature']:
			return rec
	else:
		return rec

items = [('chr2', 400, 500, 0.00250496203946)]

tItems = [t for t in items if t[3] >= thresh]

sItems = sorted(tItems, key=itemgetter(0,1,2))

x = pybedtools.BedTool(sItems)

anomolies = x.merge(d=50, c='4,4,4', o='max,min,mean')

genes = pybedtools.BedTool(infile).each(gffFilter, feature='gene')

genes.window(b=anomolies, w=1000, u=True)





