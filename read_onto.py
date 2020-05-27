# -*- coding: utf-8 -*-
"""
Created on Tue May 26 21:05:03 2020

@author: Freedom
"""

import json
import os

os.chdir("D:/ziv/cellar")

f=open("ontology/cl-simple.json","rb")
onto_data=json.load(f)
print(onto_data['graphs'][0].keys())
cells=onto_data['graphs'][0]['nodes']    # dictionary of length 2632
ids=[]
tissue=[]
cell_type=[]
tissues=['kidney', 'thymus', 'spleen', "Liver",'lymph', 'stomach', 'heart', 
         'small intestine', 'large intestine','blood',"brain","thyroid",
         'placenta','eye','other']
list=[    "Embryo", "Skeletal muscle"]
for i in cells:
    if ('lbl' in i.keys()):
        ids.append(i['id'].split('/')[-1])
        cell_type.append(i['lbl'])
        if ("kidney" in i['lbl']):
            tissue.append('kidney')
        elif ('thymocyte' in i['lbl']):
            tissue.append("thymus")
        elif ('splenic' in i['lbl']):
            tissue.append("spleen")
        elif ('hepat' in i['lbl'] or "liver" in i['lbl']):
            tissue.append("liver")
        elif ('T' in i['lbl'] or 'B' in i['lbl']):
            tissue.append('lymph')
        elif ('gastric' in i['lbl'] or 'stomach' in i['lbl']) :
            tissue.append("stomach")
        elif ('heart' in i['lbl'] or 'cardic' in i['lbl']) :
            tissue.append("heart")
        elif ('small intestine' in i['lbl']):
            tissue.append('small intestine')
        elif ('large intestine' in i['lbl']):
            tissue.append('large intestine')
        elif ('blood' in i['lbl'] or 'hema' in i['lbl']):
            tissue.append('blood')
        elif ('brain' in i['lbl']):
            tissue.append("brain")
        elif ('thyroid' in i['lbl']):
            tissue.append("thyroid")
        elif ('placenta' in i['lbl']):
            tissue.append("placenta")
        elif (' eye' in i['lbl'] or "retina" in i['lbl']):
            tissue.append("eye")
        elif ('embryo' in i['lbl']):
            tissue.append("embryo")
        elif ('muscle' in i['lbl'] or 'myo' in i['lbl']):
            tissue.append("muscle")
        else:
            tissue.append("other")    

for i in range(len(cell_type)):
    cell_type[i]=ids[i]+": "+cell_type[i]
dic={}
for i in tissue:
    dic[i]=[]
for i in range(len(cell_type)):
    dic[tissue[i]].append(cell_type[i])
