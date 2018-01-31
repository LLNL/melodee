#!/usr/bin/env python

from melodee.parser import MelodeeParser

class GeneratorFactory:
    def __init__(self):
        self.tagLookup = {}
    def register(self, function, tags):
        self.tagLookup[frozenset(tags)] = function
    def lookup(self,tags):
        return self.tagLookup[frozenset(tags)]

def main():
    import sys
    import argparse
    import os
    import glob
    import importlib
    import melodee.cardioidGenerator

    genFactory = GeneratorFactory()
    #search all the Generators and register them.
    modules = glob.glob(os.path.join(os.path.dirname(__file__),"*Generator.py"))
    for modulePath in modules:
        if os.path.samefile(modulePath,__file__):
            continue
        moduleName = os.path.split(modulePath)[-1][:-3]
        importlib.import_module("melodee."+moduleName)
        #melodee = __import__("melodee."+moduleName, globals(), locals(), [], -1)
        module = eval("melodee."+moduleName)
        for (tags, function) in module.generators.items():
            genFactory.register(function, tags)
    
    ap=argparse.ArgumentParser(description="Convert Melodee into source code")
    ap.add_argument("-g", "--generate",
                    help="Specify the type of code to generate.  Multiple options lets you specialize within a translaton method",
                    dest="generators",
                    action="append",
                    required=True,
                    )
    ap.add_argument("-o", "--output",
                    help="Change the default name of the output model.",
                    dest="output",
                    action="store",
                    )
    ap.add_argument("-t", "--target",
                    help="Specify a target model.",
                    dest="targets",
                    action="append",
                    required=True,
                    )
    ap.add_argument("filenames", nargs=argparse.REMAINDER,
                    help="Filenames containing all the targets",
                    )
    options=ap.parse_args()
    targets = options.targets
    filenames = options.filenames
    outputName = options.output
    generators = options.generators
    if not outputName:
        outputName = "_".join(targets)
    
    p = MelodeeParser()

    for filename in filenames:
        p.parse(open(filename,"r").read())
    
    if len(targets) > 1:
        targetModel = ""
        targetModel += "subsystem %s {\n" % outputName
        for model in targets:
            targetModel += "   use %s;\n" % model
        targetModel += "}\n"
        p.parse(targetModel)
        finalTarget = outputName
    else:
        finalTarget = targets[0]

    genFactory.lookup(generators)(p.getModel(finalTarget),outputName)
