#!/usr/bin/env python
"""
source : 
# https://github.com/udacity/ud120-projects/issues/331
    
convert dos linefeeds (crlf) to unix (lf)
usage: dos2unix.py
"""
import sys
original = r"C:\Users\Mehdi\GitHub\TP_Assemblage\tests\kmer.pck"
destination = r"C:\Users\Mehdi\GitHub\TP_Assemblage\tests\kmer_new.pck"
content = ''
outsize = 0
with open(original, 'rb') as infile:
    content = infile.read()
with open(destination, 'wb') as output:
    for line in content.splitlines():
        outsize += len(line) + 1
        output.write(line + str.encode('\n'))

print("Done. Saved %s bytes." % (len(content)-outsize))
