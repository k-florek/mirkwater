#!/usr/bin/env python3

import sys, os
import vcf
import json
import glob

constellation_def_path = sys.argv[1]
vcf_file_path = sys.argv[2]

threeTOone = {
    'Gly':'G','Ala':'A','Leu':'L','Met':'M','Phe':'F',
    'Trp':'W','Lys':'K','Gln':'Q','Glu':'E','Ser':'S',
    'Pro':'P','Val':'V','Ile':'I','Cys':'C','Tyr':'Y',
    'His':'H','Arg':'R','Asn':'N','Asp':'D','Thr':'T'}

#TODO:add synonmous mutations
#TODO:fix issues with +E484K sites

with open(vcf_file_path,'r') as vcffile:
    reader = vcf.Reader(vcffile)
    mutations = {}
    for record in reader:
        #Convert aa site from 3 to 1 abbv and pull out AA site
        aachange = record.INFO['ANN'][0].split('|')[10].split('.')[1]
        aaFrom = aachange[:3]
        aaTo = aachange[-3:]
        aaFrom = threeTOone[aaFrom]
        aaTo = threeTOone[aaTo]
        aaPos = aachange[3:-3]
        #add mutation to record
        mutations[record.POS] = {'ref':record.REF,'alt':record.ALT,'region':record.INFO['ANN'][0].split('|')[3],'aa_change':[aaFrom,aaPos,aaTo]}

mutation_set = []
for key in mutations.keys():
    aaformat = f"{mutations[key]['region']}:{mutations[key]['aa_change'][0]}{mutations[key]['aa_change'][1]}{mutations[key]['aa_change'][2]}"
    mutation_set.append(aaformat)

#read in all constellation defs
constellations = {}
for file in glob.glob(os.path.join(constellation_def_path,"*.json")):
    with open(file,'r') as inJSON:
        jsondata = json.load(inJSON)
        try:
            constellations[jsondata['label']] = jsondata['sites']
        except KeyError:
            constellations[jsondata['name']] = jsondata['sites']

score_table = {}
for variant in constellations.keys():
    for star in constellations[variant]:
        for mutation in mutation_set:
            if mutation.upper() == star.upper():
                try:
                    score_table[variant] += 1
                except KeyError:
                    score_table[variant] = 1

highScore = 0
result = ""
for key in score_table.keys():
    if score_table[key] > highScore:
        highScore = score_table[key]
        result = key

print(f'The best match is: {result} with {highScore} matching mutations.')
print('With the following mutations:')
print(" ".join(mutation_set))
