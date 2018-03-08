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


from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order
from melodee.cardioidGenerator import CPrintVisitor

def generateCoin(model, targetName):

    template = {'target' : targetName}
    
    givenvars = model.inputs()
    targetvars = model.varsWithAttribute("target")
    diffvars = model.diffvars()
    unknownvars = model.varsWithAttribute("unknown")
    constraintvars = model.varsWithAttribute("constraint")
    markovs = model.varsWithAttribute("markov") & diffvars
    gates = model.varsWithAttribute("gate") & diffvars
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}

    lowerVals = {var : model.info("lb",var)
                 for var in model.varsWithAttribute("lb")}
    upperVals = {var : model.info("ub",var)
                 for var in model.varsWithAttribute("ub")}
    
    good = set()        
    if model.time:
        good.add(time)

    differ = Differentiator(model, good|givenvars)
    gateJacobians = {}
    for gate in order(gates):
        (dontcare,gateJacobians[gate]) = differ.diff(diffvarUpdate[gate],gate)
    markovJacobians = { var : {} for var in markovs}
    for imarkov in order(markovs):
        for jmarkov in order(markovs):
            (dontcare,markovJacobians[imarkov][jmarkov]) = differ.diff(diffvarUpdate[imarkov],jmarkov)
    differ.augmentInstructions()
    
    dt = model.addSymbol("_dt")
    gateTargets = {}
    for gate in order(gates):
        F = model.ssa[diffvarUpdate[gate]].sympy
        L = gateJacobians[gate].sympy
        M = (F-L*gate).simplify()
        
        RLA = model.addInstruction("_%s_RLA" % gate, sympy.exp(dt*L))
        RLB = model.addInstruction("_%s_RLB" % gate, M/L*(RLA-1))
        gateTargets[gate] = (RLA,RLB)

    markovOld = {}
    for markov in order(markovs):
        markovOld[markov] = model.addSymbol("_mi_old_%s" % markov) 
    markovTargets = {}
    for imarkov in order(markovs):
        summation = 0
        for jmarkov in order(markovs):
            if imarkov == jmarkov:
                continue
            if jmarkov in markovTargets:
                thisSym = markovTargets[jmarkov]
            else:
                thisSym = markovOld[jmarkov]
            summation += markovJacobians[imarkov][jmarkov].sympy*thisSym
        sss = (diffvarUpdate[imarkov]+dt*summation)/(1-dt*markovJacobians[imarkov][imarkov].sympy)
        markovTargets[imarkov] = model.addInstruction("_mi_new_%s" % imarkov,sss)
    
    out = utility.Indenter(open(targetName+".hpp","w"))
    cprinter = CPrintVisitor(out, model.ssa, set())

    out('''
#ifndef __%(target)s_HPP__
#define __%(target)s_HPP__

#include <vector>
#include <functional>
#include "IpTNLP.hpp"

namespace %(target)s {

class ThisTNLP : public Ipopt::TNLP
{
public:
   double _dtRequested;
   double _time_lb;
   double _time_ub;''', template)
    
    out.inc()
    for var in order(givenvars)+order(targetvars):
        out('std::function<double(const double)> _func_%s;',var)
    out.dec()
    out(r'''
   
   /** default constructor */
   ThisTNLP(){}

   /** default destructor */
   virtual ~ThisTNLP(){}

   virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                             Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style);

   virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

   virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                   bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                   Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);

   virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, 
                       bool new_x, Ipopt::Number& obj_value);

   virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f){return false;}

   virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g){return false;}

   /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
   virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                           Ipopt::Number* values){return false;}

   /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
   virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                       Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                       bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                       Ipopt::Index* jCol, Ipopt::Number* values){return false;}

   virtual void finalize_solution(Ipopt::SolverReturn status,
                                  Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                                  Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                                  Ipopt::Number obj_value,
                                  const Ipopt::IpoptData* ip_data,
                                  Ipopt::IpoptCalculatedQuantities* ip_cq){}

};

void getRequirements(std::vector<std::string>&);
int getHandle(const std::string&);

Ipopt::TNLP* factory(const std::vector<std::function<double(const double)> >& functions, const double time_lb, const double time_ub, const double dt);
}

#endif''', template)

    out = utility.Indenter(open(targetName+".cpp","w"))
    cprinter = CPrintVisitor(out, model.ssa, set())

    out(r'''
#include "%(target)s.hpp"
#include <cassert>
#include <cmath>

using namespace Ipopt;
using namespace std;

namespace %(target)s {

enum ReqEnum {
   ''', template)
    out.inc()
    for var in order(givenvars)+order(targetvars):
        out('_req_%s,',var)
    out.dec()
    out(r'''
  NUMREQS
};

void getRequirements(vector<string>& reqs)
{
   reqs.clear();''', template)
    out.inc()
    for var in order(givenvars)+order(targetvars):
        out('reqs.push_back("%s");',var)
    out.dec()
    out(r'''
}

int getHandle(const string& inname)
{
   if (0) {}''', template)
    out.inc()
    for var in order(givenvars)+order(targetvars):
        out('else if (inname == "%s") { return _req_%s; }',var,var)
    out.dec()
    out(r'''

   return -1;
}

Ipopt::TNLP* factory(const vector<function<double(const double)> >& functions, const double time_lb, const double time_ub, const double dt)
{
   assert(functions.size() == NUMREQS);
   ThisTNLP* thisTNLP = new ThisTNLP();''', template)
    out.inc()
    for var in order(givenvars)+order(targetvars):
        out('thisTNLP->_func_%s = functions[_req_%s];',var,var)
    out.dec()
    out(r'''
   thisTNLP->_dtRequested = dt;
   thisTNLP->_time_lb = time_lb;
   thisTNLP->_time_ub = time_ub;
   return thisTNLP;
}

enum UnkEnum {''', template)
    out.inc()
    for var in order(unknownvars):
        out('_unknown_%s,',var)
    out.dec()
    out(r'''
   NUMUNKNOWN
};

enum ConstraintEnum {''', template)
    out.inc()
    for var in order(constraintvars):
        out('_constraint_%s,',var)
    out.dec()
    out(r'''
   NUMCONSTRAINT
};

bool ThisTNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style)
{
   n = NUMUNKNOWN;
   m = NUMCONSTRAINT;
   nnz_jac_g = n*m;
   nnz_h_lag = n*n;
   index_style = FORTRAN_STYLE;
   return true;
}

bool ThisTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u)
{''', template)
    out.inc()

    for var in order(unknownvars):
        #lb = lowerVals.get(var,'-nlp_lower_bound_inf')
        #ub = upperVals.get(var,'nlp_upper_bound_inf')
        lb = lowerVals.get(var,'-1e35')
        ub = upperVals.get(var,'1e35')
        out("x_l[_unknown_%s] = %s;",var,lb);
        out("x_u[_unknown_%s] = %s;",var,ub);

    for var in order(constraintvars):
        lb = lowerVals.get(var,'-1e35')
        ub = upperVals.get(var,'1e35')
        out("g_l[_constraint_%s] = %s;",var,lb);
        out("g_u[_constraint_%s] = %s;",var,ub);
        
    out.dec()
    out(r'''
   return true;
}

bool ThisTNLP::get_starting_point(Index , bool _init_x, Number* _x,
                                  bool _init_z, Number* _z_L, Number* _z_U,
                                  Index , bool _init_lambda, Number* _lambda)
{
   double _time = _time_lb;

   if (_init_x)
   {
''', template)
    out.inc(2)
    for var in order(givenvars):
        out('double %s = _func_%s(_time);',var,var)
    
    model.printTarget(givenvars,unknownvars,cprinter)
    
    for var in order(unknownvars):
        out('_x[_unknown_%s ] = %s;',var,var)
    out.dec(2)
    out(r'''
   }
   if (_init_z)
   {
      assert(0 && "No idea what to put for z bounds.");
   }
   if (_init_lambda)
   {
      assert(0 && "No idea how to initialize lambda.");
   }
   return true;
}

bool ThisTNLP::eval_f(Index, const Number* _x, 
                      bool _new_x, Number& _obj_value)
{
   int _tsCount = round((_time_ub-_time_lb)/_dtRequested);
   double _dt = (_time_ub-_time_lb)/_tsCount;
   _obj_value = 0;

''', template)
    out.inc()
    for var in order(unknownvars):
        out('double %s = _x[_unknown_%s];',var,var)
    model.printTarget(unknownvars,
                      diffvars,
                      cprinter)
    out.dec()
    out(r'''
   for (int _ii=0; _ii<_tsCount; _ii++)
   {
      double _time = _time_lb+_ii*_dt;''', template)
    out.inc(2)
    for var in order(givenvars):
        out('double %s = _func_%s(_time);',var,var)

    out("//get the targets")
    good |= unknownvars|givenvars|diffvars
    targetSet = model.allDependencies(good,targetvars)-good
    model.printSet(targetSet,cprinter)
    good |= targetSet


    out("//get the gate updates (diagonalized exponential integrator)")
    gateGoals = set()
    for RLA,RLB in gateTargets.values():
        gateGoals.add(RLA)
        gateGoals.add(RLB)
    good.add(dt)
    gateSet = model.allDependencies(good,gateGoals)-good
    model.printSet(gateSet,cprinter)
    good |= gateSet

    out("//get the other differential updates")
    diffGoals = set([diffvarUpdate[var] for var in diffvars-gates])
    diffSet = model.allDependencies(good,diffGoals)-good
    model.printSet(diffSet, cprinter)
    good |= diffSet
    
    out("//Do the markov update (1 step rosenbrock with gauss siedel)")
    for var in order(markovs):
        out("double %s = %s;",markovTargets[var],diffvarUpdate[var])
    out("int _count=0;")
    out("double _error;")
    out("do")
    out("{")
    out.inc()
    cprinter.pushStack()


    for var in order(markovs):
        out("double %s = %s;",markovOld[var], markovTargets[var])
    markovGoals = set(markovOld.values())
    good |= markovGoals
    cprinter.declaredStack[-1] |= set(["%s" % var for var in markovTargets.values()])
    markovSet = model.allDependencies(good,set(markovTargets.values()))
    model.printSet(markovSet,cprinter)

    out("_error = 0;")
    for var in order(markovs):
        old = markovOld[var]
        new = markovTargets[var]
        out("_error += (%s-%s)*(%s-%s);",old,new,old,new)
    out("_count++;")

    cprinter.popStack()
    out.dec()
    out("} while (_error > 1e-100 && _count<50);")
    
    for var in order(diffvars-gates-markovs):
        out("%s += _dt*%s;",var,diffvarUpdate[var])

    for var in order(gates):
        RLA,RLB = gateTargets(var)
        out("%s = %s*%s + %s;",var,RLA,var,RLB)

    for var in order(markovs):
        out("%s += _dt*%s;",var,markovTargets[var])
        
    for var in order(targetvars):
        out("double _residual_%s = %s - _func_%s(_time);",var,var,var)
        out("_obj_value += _dt*_residual_%s*_residual_%s;",var,var)
    out.dec(2)
    out(r'''
   }
   return true;
}

} // namespace HergCoinStandalone
''' % template)

generators = {
    frozenset(["coin"]) : generateCoin
}
