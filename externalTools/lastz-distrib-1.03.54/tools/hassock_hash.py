#!/usr/bin/env python
"""
Python implementation of the hash used for "hashed sequence index" files.

The "hassock" hash is a variant of Austin Appleby's MurmurHash2.  The
latter is described (as of Apr/2009) at
	murmurhash.googlepages.com
This variant is based on the endian-neutral version found at
	murmurhash.googlepages.com/MurmurHashNeutral2.cpp
and differs in the following ways:
	(a)	The "seed" is hardwired.
	(b)	We parse the data block in reverse;  this allows the caller to
		prepend an additional seed pattern to his buffer, potentially
		getting better mixing for the bits in the final incorporated
		bytes.
	(c)	The last three bytes are incorporated in a different order than
		they were in MurmurHash2, because the code just works out better
		this way.

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys

seed = 0x5C3FC4D3
mult = 0x87C10417

def hassock_hash(s):
	ix = len(s)
	h  = seed ^ ix							# h = seed ^ len;
	while (ix >= 4):
		k =  ord(s[ix-1])					# k  = *(--data);
		k |= ord(s[ix-2]) << 8				# k |= *(--data) << 8;
		k |= ord(s[ix-3]) << 16				# k |= *(--data) << 16;
		k |= ord(s[ix-4]) << 24				# k |= *(--data) << 24;

		k = (k * mult) & 0xFFFFFFFF			# k *= m; 
		k ^= k >> 24 						# k ^= k >> r;
		k = (k * mult) & 0xFFFFFFFF			# k *= m; 

		h = (h * mult) & 0xFFFFFFFF			# h *= m; 
		h ^= k								# h ^= k;
		ix -= 4

	if (ix >= 3):
		h ^= ord(s[2]) << 16				# h ^= *(--data) << 16;
	if (ix >= 2):
		h ^= ord(s[1]) << 8					# h ^= *(--data) << 8;
	if (ix >= 1):
		h ^= ord(s[0])						# h ^= *(--data);
		h = (h * mult) & 0xFFFFFFFF			# h *= m; 

	h ^= h >> 13							# h ^= h >> 13;
	h = (h * mult) & 0xFFFFFFFF				# h *= m; 
	h ^= h >> 15							# h ^= h >> 15;

	return h


# main program to test

def main():
	m = None

	strings = []

	for s in sys.argv[1:]:
		if (s.startswith("--mod=")): m = int(s.split("=",1)[1])
		else:                        strings += [s]

	if (strings != []):
		for s in strings:
			demonstrate_hash(s,m)
	else:
		for line in sys.stdin:
			line = line.rstrip()
			demonstrate_hash(line,m)

def demonstrate_hash(s,m):
	if (m == None): print "%08X: %s" % (hassock_hash(s),s)
	else:           print "%d: %s"   % (hassock_hash(s)%m,s)


if __name__ == "__main__": main()

