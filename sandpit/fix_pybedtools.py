#fix pybedtools fuckery

import pybedtools
from operator import itemgetter
import re

windows = [('chr_1',100,200,5.0),('chr_10',150,250,5.0),('chr_2',200,350,1.0)]

stringyWindows = [ (w,x,y,str(z)) for w,x,y,z in windows ]

anomaliesBED = pybedtools.BedTool(stringyWindows)


def natural_sort(l,key=str): 
    """Sorting of weird scaffold names for chromosome painting."""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda item: [ convert(c) for c in re.split('([0-9]+)', key(item)) ] 
    return sorted(l, key = alphanum_key)


print natural_sort(anomaliesBED, key=itemgetter(0))


###
a = list()
b = [1,2]

if len(a) > 0 and len(b) > 0:
	print 'passed,yo!'