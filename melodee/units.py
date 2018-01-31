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

import math

class Unit:
    def __init__(self, unitSystem, bases, scale):
        self.system = unitSystem
        self.bases = bases
        self.scale = scale

    def isCompatibleWith(self, other):
        if other == None:
            return False
        if self.system != other.system:
            return False
        for (base,power) in self.bases.items():
            if power != 0:
                if base not in other.bases:
                    return False
                elif other.bases[base] != power:
                    return False
        for (base,power) in other.bases.items():
            if power != 0:
                if base not in self.bases:
                    return False
        return True
    def __eq__(self, other):
        return self.isCompatibleWith(other) and self.scale == other.scale
    def __ne__(self, other):
        return not self==other
    def copy(self):
        basesCopy = dict(self.bases)
        return Unit(self.system, basesCopy, self.scale)
    def __mul__(self, other):
        assert(self.system == other.system)
        ret = self.copy()
        for (base,power) in other.bases.items():
            ret.bases[base] = ret.bases.get(base, 0) + power
        ret.scale += other.scale
        return ret
    def __div__(self, other):
        return self.__truediv__(other)
    def __truediv__(self, other):
        assert(self.system == other.system)
        ret = self.copy()
        for (base,power) in other.bases.items():
            ret.bases[base] = ret.bases.get(base, 0) - power
        ret.scale -= other.scale
        return ret
    def __pow__(self, number):
        ret = self.copy()
        for base in self.bases.keys():
            ret.bases[base] *= number
        ret.scale = ret.scale*number
        return ret
    def __str__(self):
        return self.system.strify(self)
    def __repr__(self):
        return 'Unit(si,'+repr(self.bases)+','+repr(self.scale)+')'
    def convertTo(self, other):
        '''
        Gives factor so that self * factor = providedUnit
        Test mentally with millimeters and meters
        '''
        assert self.isCompatibleWith(other)
        return 10**(self.scale-other.scale)

class UnitSystem:
    def __init__(self, possibleBases):
        self.knownUnits = {}
        self.possibleBases = set()
        for base in possibleBases:
            self.addBase(base)
    def addBase(self, base):
        if base not in self.possibleBases:
            self.possibleBases.add(base)
            self.register(base, Unit(self, {base : 1}, 0))
    def register(self, name, unit):
        self.knownUnits[name] = unit
    def scaleToUnit(self,scale):
        return Unit(self, {}, math.log10(scale))
    def get(self, unitName):
        return self.knownUnits[unitName]
    def strify(self, unit):
        '''find the smallest unit expr by greedy gradient descent'''
        def unitDot(xxx,yyy):
            allBases = set(xxx.bases.keys()) | set(yyy.bases.keys())
            dotp = xxx.scale*yyy.scale
            for base in allBases:
                dotp += xxx.bases.get(base, 0)*yyy.bases.get(base, 0)
            return dotp
        numerator = []
        denominator = []
        remainder = unit
        while unitDot(remainder,remainder) > 0:
            maxUnitName = None
            maxUnitDot = 0
            for (thisName,thisUnit) in self.knownUnits.items():
                thisSize = math.sqrt(unitDot(thisUnit,thisUnit))
                if thisSize==0:
                    continue
                thisUnitDot = unitDot(remainder,thisUnit)/thisSize
                if (maxUnitName == None or
                    abs(thisUnitDot) > abs(maxUnitDot) or
                    (abs(thisUnitDot) == abs(maxUnitDot) and
                     len(thisName) < len(maxUnitName)
                     )):
                    maxUnitDot = thisUnitDot
                    maxUnitName = thisName
            assert(maxUnitName != None)
            if maxUnitDot > 0:
                numerator.append(maxUnitName)
                remainder = remainder / self.get(maxUnitName)
            else:
                denominator.append(maxUnitName)
                remainder = remainder * self.get(maxUnitName)
        if not numerator and not denominator:
            return "1"
        elif numerator and not denominator:
            return "*".join(numerator)
        elif not numerator and denominator:
            return "1/"+"/".join(denominator)
        else:
            return "*".join(numerator)+"/"+"/".join(denominator)

    def __getattr__(self, name):
        try:
            return self.get(name)
        except KeyError:
            raise AttributeError(name)
    
class Si(UnitSystem):
    def __init__(self):
        baseSi = ['mol', 'm', 's', 'A', 'g', 'cd', 'K', 'bool']
        UnitSystem.__init__(self, baseSi)

        s = list(baseSi)
        l = []

        x='unitless'; self.register(x, self.scaleToUnit(1))
        x='dimensionless'; self.register(x, self.unitless)
        x='1'; self.register(x, self.unitless)

        x='L'; s.append(x); self.register(x, self.m**3 * self.scaleToUnit(1e-3))
        x='liter'; l.append(x); self.register(x, self.L)
        x='litre'; l.append(x); self.register(x, self.L)
        x='M'; s.append(x); self.register(x, self.mol / self.L)
        x='molar'; l.append(x); self.register(x, self.M)
        x='C'; s.append(x); self.register(x, self.A*self.s)
        x="coulumb"; l.append(x); self.register(x, self.C)
        x='V'; s.append(x); self.register(x, self.scaleToUnit(1e3)*self.m**2*self.g/(self.s**3*self.A))
        x='volt'; l.append(x); self.register(x, self.V)
        x='ohm'; l.append(x); self.register(x, self.V/self.A)
        x='S'; s.append(x); self.register(x, self.unitless/self.ohm)
        x='siemen'; l.append(x); self.register(x, self.S)
        x='F'; s.append(x); self.register(x, self.C/self.V)
        x='farad'; l.append(x); self.register(x, self.F)
        x='N'; s.append(x); self.register(x, self.scaleToUnit(1e3)*self.m*self.g/self.s**2)
        x='newton'; l.append(x); self.register(x, self.N)
        x='J'; s.append(x); self.register(x, self.N*self.m)
        x='joule'; l.append(x); self.register(x,self.J)
        x='Pa'; s.append(x); self.register(x, self.N/self.m**2)
        x='pascal'; l.append(x); self.register(x,self.Pa)
        x='W'; s.append(x); self.register(x, self.J/self.s)
        x='watt'; l.append(x); self.register(x, self.W)
        x='Hz'; s.append(x); self.register(x, self.unitless/self.s)
        x='hertz'; s.append(x); self.register(x, self.Hz)
        
        x='meter'; l.append(x); self.register(x, self.m)
        x='metre'; l.append(x); self.register(x, self.m)
        x='second'; l.append(x); self.register(x, self.s)
        x='ampere'; l.append(x); self.register(x, self.A)
        x='amp'; l.append(x); self.register(x, self.A)
        x='gram'; l.append(x); self.register(x, self.g)
        x='candela'; l.append(x); self.register(x, self.cd)
        x='kelvin'; l.append(x); self.register(x, self.K)

        x='eV'; s.append(x); self.register(x, self.J*self.scaleToUnit(6.241509e-18))
        x='electron_volt'; l.append(x); self.register(x, self.eV)
        x='cal'; s.append(x); self.register(x,self.J*self.scaleToUnit(0.239006))
        x='calorie'; l.append(x); self.register(x,self.cal)
        x='Ang'; s.append(x); self.register(x, self.m*self.scaleToUnit(1e-10))
        x='angstrom'; l.append(x); self.register(x, self.Ang)
        
        shortNames = s
        longNames = l

        shortPrefixes = {
            'y': 1e-24,
            'z' : 1e-21,
            'a' : 1e-18,
            'f' : 1e-15,
            'p' : 1e-12,
            'n' : 1e-9,
            'u' : 1e-6,
            'm' : 1e-3,
            'c' : 1e-2,
            'd' : 1e-1,
            'da' : 1e1,
            'h' : 1e2,
            'k' : 1e3,
            'M' : 1e6,
            'G' : 1e9,
            'T' : 1e12,
            'P' : 1e15,
            'E' : 1e18,
            'Z' : 1e21,
            'Y' : 1e24,
        }
        for shortName in shortNames:
            for (prefix,scale) in shortPrefixes.items():
                self.register(prefix+shortName, self.get(shortName)*self.scaleToUnit(scale))

        longPrefixes = {
            'yotto': 1e-24,
            'zepto' : 1e-21,
            'atto' : 1e-18,
            'femto' : 1e-15,
            'pico' : 1e-12,
            'nano' : 1e-9,
            'micro' : 1e-6,
            'milli' : 1e-3,
            'centi' : 1e-2,
            'deci' : 1e-1,
            'deca' : 1e1,
            'hecto' : 1e2,
            'kilo' : 1e3,
            'mega' : 1e6,
            'giga' : 1e9,
            'tera' : 1e12,
            'peta' : 1e15,
            'exa' : 1e18,
            'zetta' : 1e21,
            'yotta' : 1e24,
        }
        for longName in longNames:
            for (prefix, scale) in longPrefixes.items():
                self.register(prefix+longName, self.get(longName)*self.scaleToUnit(scale))

        self.register('cc', self.cm**3)

if __name__=='__main__':
    si  = Si()
    assert (si.mV/si.ms).isCompatibleWith(si.uA/si.uF)
    print((si.mV/si.ms).convertTo(si.uA/si.uF)-1)
    assert (si.mV/si.ms).convertTo(si.uA/si.uF) == 1
    
