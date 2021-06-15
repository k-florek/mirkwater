#!/usr/bin/env python3

import sys, os
import vcf
import json
import glob
import csv

constellation_def_path = sys.argv[1]
vcf_file_path = sys.argv[2]

threeTOone = {
    'Gly':'G','Ala':'A','Leu':'L','Met':'M','Phe':'F',
    'Trp':'W','Lys':'K','Gln':'Q','Glu':'E','Ser':'S',
    'Pro':'P','Val':'V','Ile':'I','Cys':'C','Tyr':'Y',
    'His':'H','Arg':'R','Asn':'N','Asp':'D','Thr':'T'}

#TODO:add synonmous mutations

with open(vcf_file_path,'r') as vcffile:
    reader = vcf.Reader(vcffile)
    mutations = {}
    sampleName = reader.samples[0]
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
class Constellation:
    def __init__(self,label,sites,alt_rule=None):
        self.label = label
        self.sites = sites
        self.alt_rule = alt_rule

constellations = []
for file in glob.glob(os.path.join(constellation_def_path,"*.json")):
    with open(file,'r') as inJSON:
        jsondata = json.load(inJSON)
        alt_rule = None
        for rules in jsondata['rules']:
            if jsondata['rules'][rules] == 'alt':
                alt_rule = rules
                constellations.append(Constellation(jsondata['label'], jsondata['sites'],alt_rule))
        if alt_rule == None:
            constellations.append(Constellation(jsondata['label'], jsondata['sites']))

class Score:
    def __init__(self,constellation,alt_rule = None):
        self.label = constellation
        self.score = 0
        self.mutations = []
        self.alt_rule = alt_rule
        self.alt_rule_hit = False
        if alt_rule != None:
            self.alt_rule_required = True
        else:
            self.alt_rule_required = False

    def hit(self,mutation):
        self.mutations.append(mutation)
        if self.alt_rule_required:
            if mutation == self.alt_rule:
                self.alt_rule_hit = True
        self.score += 1


scores = []
for constellation in constellations:
    score = Score(constellation.label,constellation.alt_rule)
    for site in constellation.sites:
        for mutation in mutation_set:
            if mutation.upper() == site.upper():
                score.hit(mutation.upper())
    scores.append(score)

highScore = None
highScorePoints = 0
for score in scores:
    if score.alt_rule_required and score.alt_rule_hit:
        if score.score > highScorePoints:
            highScore = score
            highScorePoints = score.score
    if not score.alt_rule_required:
        if score.score > highScorePoints:
            highScore = score
            highScorePoints = score.score

mutationHitList = '; '.join(highScore.mutations)
with open(f'{sampleName}.results.csv','w') as outfile:
    outfile.write(f'{sampleName},{highScore.label},{mutationHitList}\n')

print('---Hit List---')
for score in scores:
    print('--------------')
    print(score.label)
    print(score.mutations)
    print(score.score)

print('\n---Best Hit---')
print(highScore.label)
print(highScore.mutations)
print(highScorePoints)
