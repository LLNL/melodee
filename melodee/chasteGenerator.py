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
from sympy.printing.c import C99CodePrinter
from sympy.core import S


from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order

from melodee.cardioidGenerator import repeat,MyCCodeSympyPrinter,CPrintVisitor

def pretty(symbol):
    return str(symbol)

def generateChaste(model, targetName):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    params = model.varsWithAttribute("param")
    gates = model.varsWithAttribute("gate") & diffvars
    tracevars = model.varsWithAttribute("trace")
    markovs = model.varsWithAttribute("markov") & diffvars

    lbLookup = {}
    ubLookup = {}
    for var in gates|markovs:
        lbLookup[var] = 0
        ubLookup[var] = 1

    for var in model.varsWithAttribute("lb"):
        lbLookup[var] = model.info("lb",var)
        ubLookup[var] = model.info("ub",var)
    
    V = model.input("V")
    V_init = model.output("V_init")
    Iion = model.output("Iion")
    hasCai = False
    if "Cai" in model.outputs():
        Cai = model.output("Cai")
        hasCai = True
    
    out = utility.Indenter(open(targetName+".hpp","w"))
    out('''
#ifndef %(target)s_HPP_
#define %(target)s_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"
#include <vector>

/**
 * This class represents the Luo-Rudy 1991 system of equations,
 * with support for being compiled into a .so and loaded at run-time.
 */
class %(target)s : public AbstractCardiacCell, public AbstractDynamicallyLoadableEntity
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
        archive & boost::serialization::base_object<AbstractDynamicallyLoadableEntity>(*this);
    }

    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();

public:
    /**
     * Constructor
     *
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    %(target)s(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
               boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~%(target)s();

    /**
     * Fill in a vector representing the RHS of the Luo-Rudy 1991 system
     * of Odes at each time step, y' = [y1' ... yn'].
     * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Returns the ionic current
     *
     * @param pStateVariables  optional state at which to evaluate the current
     * @return the total ionic current
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);''',template)
    if hasCai:
        out('''
    /**
     * Get the intracellular calcium concentration
     *
     * @return the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration();
''')
    out(r'''
};

#include "SerializationExportWrapper.hpp"
    CHASTE_CLASS_EXPORT(%(target)s)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a %(target)s instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const %(target)s * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
    const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
    ar << p_solver;
    ar << p_stimulus;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a %(target)s instance (using existing constructor).
 *
 * NB this constructor allocates memory for the other member variables too.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, %(target)s * t, const unsigned int file_version)
{

    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    ar >> p_solver;
    ar >> p_stimulus;
    ::new(t)%(target)s(p_solver, p_stimulus);
}
}
} // namespace ...

#endif // %(target)s_

''', template)

    out = utility.Indenter(open(targetName+".cpp","w"))
    out(r'''
#include "%(target)s.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>
//#include <iostream>
#include "Exception.hpp"

enum StateVarEnum {
''',template)
    for var in order(diffvars):
        out("_enum_%s,",var)
    out(r'''
  NUMSTATES
};

/**
 * Constructor
 */
%(target)s::%(target)s(
    boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
    boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(pSolver, NUMSTATES+1, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<%(target)s>::Instance();

    Init();
}

/**
 * Destructor
 */
%(target)s::~%(target)s(void)
{
}
''',template)
    if hasCai:
        out('''
double %s::GetIntracellularCalciumConcentration()
{
    return mStateVariables[_enum_%s+1];
}
''', targetName, Cai)
    out(r'''

/**
 * Fill in a vector representing the RHS of the %(target)s system
 * of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param time  the current time, in milliseconds
 * @param rY  current values of the state variables
 * @param rDY  to be filled in with derivatives
 */
void %(target)s::EvaluateYDerivatives(double time,
                                                      const std::vector<double> &_rY,
                                                      std::vector<double> &_rDY)
{
''', template)
    out.inc()

    out("double %s = _rY[0];", V)
    for var in order(diffvars):
        out("double %s = _rY[_enum_%s+1];", var, var)
    ii=0
    for var in order(params):
        out("double %s = GetParameter(%d);", var, ii);
        ii+=1
        
    cprinter = CPrintVisitor(out, model.ssa, params)
    good=diffvars|set([V])|params
    model.printTarget(good, set(diffvarUpdate.values())|set([Iion]), cprinter)

    out("_rDY[0] = -%s;", Iion)
    for var in order(diffvars):
        out("_rDY[_enum_%s+1] = %s;", var, diffvarUpdate[var])
    out.dec()
    out('''
}


double %(target)s::GetIIonic(const std::vector<double>* pStateVariables)
{
''', template)
    out.inc()

    out("double %s = _rY[0];", V)
    for var in order(diffvars):
        out("double %s = _rY[_enum_%s+1];", var, var)
    ii=0
    for var in order(params):
        out("double %s = GetParameter(%d);", var, ii);
        ii+=1

    cprinter = CPrintVisitor(out, model.ssa, params)
    model.printTarget(good, set([Iion]), cprinter)

    out("assert(!std::isnan(%s));",Iion)
    out("return %s;",Iion)

    out.dec()
    out('''
}

void %(target)s::VerifyStateVariables()
{
    const std::vector<double>& _rY = rGetStateVariables();
''', template)
    out.inc()

    for var in order(diffvars):
        out("double %s = _rY[_enum_%s+1];", var, var)
        if var in lbLookup:
            out("if (%s < %s)", var, lbLookup[var])
            out("{")
            out.inc()
            out('EXCEPTION(DumpState("%s has gone under its lower bound. Check model parameters, for example spatial stepsize"));', var)
            out.dec()
            out("}")
        if var in ubLookup:
            out("if (%s > %s)", var, ubLookup[var])
            out("{")
            out.inc()
            out('EXCEPTION(DumpState("%s has gone over its upper bound. Check model parameters, for example spatial stepsize"));', var)
            out.dec()
            out("}")

    out.dec()
    out('''
}

template<>
void OdeSystemInformation<%(target)s>::Initialise(void)
{
''', template)
    out.inc()

    cprinter = CPrintVisitor(out, model.ssa, params)
    model.printTarget(set(),diffvars|set([V_init])|params,cprinter)
    
    out('this->mVariableNames.push_back("%s");',V)
    out('this->mVariableUnits.push_back("mV");') # FIXME, make generic
    out('this-mInitialConditions.push_back(%s);',V_init)

    for var in order(diffvars):
        out('this->mVariableNames.push_back("%s");',var)
        out('this->mVariableUnits.push_back("%s");', model.ssa[var].astUnit.rawUnit)
        out('this-mInitialConditions.push_back(%s);',var)

    for var in order(params):
        out('this->mParameterNames.push_back("%s");',var)
        out('this->mParameterUnits.push_back("");') #FIXME, put in a unit if we have one.
        out('this->mParameters.push_back(%s);',var)
        
    out.dec()
    out('''
    this->mInitialised = true;
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(%(target)s)

extern "C" {
    /**
     * C-style constructor with a standard name, so the main program can create instances
     * of this cell model when it is loaded from a .so.
     *
     * @param pSolver  ODE solver used to simulate the cell
     * @param pIntracellularStimulus  intracellular stimulus
     */
    AbstractCardiacCellInterface* MakeCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                         boost::shared_ptr<AbstractStimulusFunction> pStimulus){
         return new %(target)s(pSolver, pStimulus);
    }
}
''', template)

generators = {
    frozenset(["chaste"]) : generateChaste
}
