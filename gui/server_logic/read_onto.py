# -*- coding: utf-8 -*-
"""
Created on Tue May 26 21:05:03 2020

@author: Freedom
"""

import json



def get_dic():
    f=open("markers/cl-simple.json","rb")
    onto_data=json.load(f)
    cells=onto_data['graphs'][0]['nodes']    # dictionary of length 2632
    tissue=[]
    cell_type=[]
    tissues=['kidney', 'thymus', 'spleen', "liver",'lymph', 'stomach', 'heart',
         'small intestine', 'large intestine','blood',"brain","thyroid",
         'placenta','eye','other', "embryo", "muscle"]

    identifiers=[['kidney'], ['thymocyte'], ['splenic'], ["hepat","liver"],
            ['T','B'], ['stomach','gastric'], ['heart','cardiac'],
         ['small intestine'], ['large intestine'],['blood'],["brain"],["thyroid"],
         ['placenta'],['eye','retina'], ['other'],["embryo"], ["muscle",'myo']]

    for i in cells:
        f=0
        if ('lbl' in i.keys()):
            for j in range(len(identifiers)):
                for k in identifiers[j]:
                    if ((k in i['lbl']) and f==0):
                        tissue.append(tissues[j])
                        f=1
            if (f==0):
                tissue.append('other')
            cell_type.append(i['id'].split('/')[-1]+': '+i['lbl'])
        else:
            tissue.append('other')
            cell_type.append(i['id'].split('/')[-1])
    dic={}
    for i in tissue:
        dic[i]=[]
    for i in range(len(cell_type)):
        dic[tissue[i]].append(cell_type[i])
    return dic
