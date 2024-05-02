#!/bin/env python3
""" generation of module dependances for a xxx.f90 or xxx.F90 file """

import string,sys
from subprocess import Popen, PIPE

def noComment(mot):
    """ remove comment characters in a word """
    ind= mot.find('!')
    if ind != -1:
        mot= mot[:ind]
    return mot

class md9:
    """ a file, its module and its dependances """
    def __init__(self,fic):
        self.name= fic
        self.used= []
        #find all module dependances
        listeLignes= open(self.name,'r',encoding='utf8',errors='ignore').readlines()
        self.module= "dummy"
        for ligne in listeLignes:
            ligne= ligne.lower()
            if ligne.count('use'):
                ligne= ligne.replace(',',' ')
                listeMots= ligne.split()
                if listeMots[0] == 'use':
                    mod= noComment(listeMots[1])
                    if mod not in self.used and mod != 'intrinsic':
                        self.used.append(mod)
            elif ligne.count('module'):
                listeMots= ligne.split()
                if listeMots[0] == 'module' and \
                   noComment(listeMots[1]) != 'procedure':
                    self.module= listeMots[1]
    def utilise(self,other):
        return (other.module in self.used)

#list the file names
listeClasse= list()
for x in sys.argv[1:]:
    if x.endswith('.f90') or x.endswith('.F90'): listeClasse.append(md9(x))
#sort the file list
listeRes= []
while listeClasse:
    lgav= len(listeClasse)
    for fic in listeClasse:
        if len(list(filter(lambda x,y=fic:y.utilise(x),listeClasse))) == 0:
            listeRes.append(fic)
            listeClasse.remove(fic)
    lgap= len(listeClasse)
    if lgav == lgap and lgap != 0:
        raise RuntimeError("make_depend: cross-references found")
#result output
resstr=''
for fic in listeRes:
    resstr=resstr+fic.name+' '
print(resstr)
