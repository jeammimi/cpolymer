# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 15:24:01 2014

@author: jarbona
"""
description_pair = {"lj/cut":{"template":"{epsilon:.2f} {sigma:.2f} {cutoff1:.2f}",
                                "required":["idpair1","idpair2","epsilon","sigma"],
                                "optional":{"cutoff1":""}},
                    "lj/expand":{"template":"{epsilon:.2f} {sigma:.2f} {delta:.2f} {cutoff1:.2f}",
                                "required":["idpair1","idpair2","epsilon","delta","sigma"],
                                "optional":{"cutoff1":""}},
                    "soft":{"template":"{A:.2f} {rc:.2f}",
                                                    "required":["idpair1","idpair2","A"],
                                                    "optional":{"rc":""}}
                    }
                    

class Pair:
    def __init__(self,typep,hybrid=False,**args):
        self.typep = typep
        self.hybrid = hybrid
        
        if typep not in list(description_pair.keys()):
            raise typep + "not described"
        self.args = args
        
        
        for k in description_pair[typep]["required"]:
            if k not in list(args.keys()):
                print(k) 
                print(" argument needed")
                raise 
        if "optional" in description_pair[typep]:
            for k,v in list(description_pair[typep]["optional"].items()):
                if k not in list(args.keys()):
                    args[k] = v
        
        
    def __repr__(self):
        hybrid = " "
        
        if self.hybrid:
            hybrid = " {0} ".format(self.typeb)
        return "pair_coeff {0} {1} {2}".format(self.args["idpair1"],self.args["idpair2"],hybrid) + description_pair[self.typep]["template"].format(**self.args)

if __name__ == "__main__":
    print((Pair("lj/cut",idpair1=1,idpair2=2,epsilon=1.3,sigma=1.1)))

    #print Bond("harmonic",idbond=1,distance=1.4,hybrid=True)
