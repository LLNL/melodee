#!/usr/bin/env python
#### <license> ####
# Copyright (c) 2016-2017, Lawrence Livermore National Security,
# LLC. Produced at the Lawrence Livermore National Laboratory. Written
# by Robert Blake <blake14@llnl.gov>.
# 
# LLNL-CODE-720003.
# All rights reserved.
#
# This file is part of MELODEE. For details, see
# http://github.com/llnl/melodee.
#
# Licensed under the Apache License, Version 2.0 (the "Licensee"); you
# may not use this file except in compliance with the License.  You may
# obtain a copy of the License at:
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing
# permissions and limitations under the license.
#### </license> #### 

import sys
import re
import sympy

class Indenter:
    def __init__(self, fff=sys.stdout, indentString="   "):
        self.indent = indentString
        self.indentAmount = 0
        self.outfile = fff

    def __call__(self, string, *args, **kwargs):
        initialstrip = string.lstrip('\n')
        if not initialstrip:
            print >>self.outfile, string,
            return
        indentString = self.indent*self.indentAmount
        for line in initialstrip.split('\n'):
            outline = indentString+line
            if line == "":
                print >>self.outfile, "\n",
            elif kwargs:
                print >>self.outfile, (outline % kwargs)
            elif len(args) == 1 and type(args[0]) == dict:
                print >>self.outfile, (outline % args[0])
            elif len(args) == 1:
                print >>self.outfile, (outline % args[0])
            elif len(args) == 0:
                print >>self.outfile, outline
            else:
                print >>self.outfile, (outline % args)

    def inc(self, indentAmount=1):
        self.indentAmount += indentAmount

    def dec(self, indentAmount=1):
        self.indentAmount -= indentAmount

def order(iterable):
    ret = [item for item in iterable]
    if ret and isinstance(ret[0], sympy.Symbol):
        ret.sort(key=str)
    else:
        ret.sort()
    return ret

def unzip(lll):
    """Does the opposite of zip"""
    return zip(*lll)
