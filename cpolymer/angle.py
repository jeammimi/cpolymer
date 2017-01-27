# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 15:24:01 2014

@author: jarbona
"""
description_angle = {"harmonic":{"template":"{K} {theta}",
                                "required":["idangle","K","theta"]},
                    "cosine/delta":{"template":"{K} {theta}",
                                "required":["idangle","K","theta"]}
                }
                    

class Angle:
    def __init__(self,typea,hybrid=False,**args):
        self.typea = typea
        self.hybrid = hybrid
        
        if typea not in list(description_angle.keys()):
            raise typea + "not described"
        self.args = args
        
        
        for k in description_angle[typea]["required"]:
            if k not in list(args.keys()):
                print(k) 
                print( " argument needed")
                raise
        if "optional" in description_angle[typea]:
            for k,v in list(description_angle[typea]["optional"].items()):
                if k not in list(args.keys()):
                    args[k] = v
        
        
    def __repr__(self):
        hybrid = " "
        
        if self.hybrid:
            hybrid = " {0} ".format(self.typeb)
        return "angle_coeff {0} {1}".format(self.args["idangle"],hybrid) + description_angle[self.typea]["template"].format(**self.args)

if __name__ == "__main__":
    print((Angle("harmonic",idangle=1,K=80,theta=180)))

    #print Bond("harmonic",idbond=1,distance=1.4,hybrid=True)
