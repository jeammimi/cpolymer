# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 15:24:01 2014

@author: jarbona
"""
description_bond = {"harmonic":{"template":"{K:.2f} {R0:.2f}",
                                "required":["idbond","K","R0"]},
                    "fene":{"template":"{K:.2f} {R0:.2f} {epsilon:.2f} {sigma:.2f}",
                                "required":["idbond","K","R0","epsilon","sigma"]},
                    "fene/expand":{"template":"{K:.2f} {R0:.2f} {epsilon:.2f} {sigma:.2f} {delta:.2f}",
                                "required":["idbond","K","R0","epsilon","sigma","delta"]}}

class Bond:
    def __init__(self,typeb,hybrid=False,**args):
        self.typeb = typeb
        if typeb not in list(description_bond.keys()):
            print((typeb + "not described"))
            raise typeb + "not described"
        self.args = args
        self.hybrid = hybrid

        
        for k in description_bond[typeb]["required"]:
            if k not in list(args.keys()):
                print(k) 
                print(" argument needed")
                raise 
        if "optional" in description_bond[typeb]:
            for k,v in list(description_bond[typeb]["optional"].items()):
                if k not in list(args.keys()):
                    args[k] = v
        
    def __repr__(self):
        hybrid = " "
        if self.hybrid:
            hybrid = " {0} ".format(self.typeb)
        return "bond_coeff {0}{1}".format(self.args["idbond"],hybrid) + description_bond[self.typeb]["template"].format(**self.args)

if __name__ == "__main__":
    print((Bond("harmonic",idbond=1,R0=1.4,K=1)))
    print((Bond("harmonic",idbond=1,R0=1.4,K=1,hybrid=True)))
    print((Bond("fene",idbond=1,R0=1.4,K=1,epsilon=1,sigma=1,hybrid=True)))
    #print Bond("harmonic",idbond=1,distance=1.4,hybrid=True)
