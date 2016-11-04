#!/usr/bin/env python2.7

import os
import argparse
import sys

def main():
    print "Got command line arguments: %s" % sys.argv
    for line in sys.stdin:
        print line

if __name__=="__main__":
    main()
