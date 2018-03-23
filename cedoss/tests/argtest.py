#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 13:37:33 2018

tests command line python implementation

@author: chris
"""

import sys

def greetings(word1=sys.argv[1], word2=sys.argv[2]):
    print("{} {}".format(word1, word2))

if __name__ == "__main__":
    greetings()
    greetings("Bonjour", "monde")