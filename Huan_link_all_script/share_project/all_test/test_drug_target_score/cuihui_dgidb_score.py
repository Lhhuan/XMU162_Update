import pandas as pd
import requests
import json
import os
import numpy as np

def SearchDrugGeneScore(drug):
    key = {"drugs":drug}
    tem = requests.get('http://dgidb.org/api/v2/interactions.json?',params=key)
    scores = tem.json()
    b = [['Drugs','interactionTypes','geneName','sources','pmids','score']]
    if scores["matchedTerms"] == []:
        b = []
    else:
        if scores["matchedTerms"][0]['searchTerm'] == drug:
            for i in range(len(scores["matchedTerms"])):
                for item in scores["matchedTerms"][i]['interactions']:
                    a = [scores["matchedTerms"][i]['searchTerm'],",".join(item['interactionTypes']),
                         item['geneName'],",".join(item['sources']),str(item['pmids']).strip('[]'),str(item['score'])]
                    b= np.vstack((b,a))
            if len(b)==1:
                b = []
        else:
            b = []
    return b[1:]

file = pd.read_table('uniq_drug.tsv',sep='\t')
drug = file.loc[:,['drug_claim_name','drug_claim_primary_name','drug_name','drug_chembl_id']]
with open('sort_unique_no_indication_drug_target_score.txt','w') as f:
    title = 'drug_claim_name\tdrug_claim_primary_name\tdrug_name\tdrug_chembl_id\tinteractionTypes\tgeneName\tsources\tpmids\tscore\ttag\n'
    f.writelines(title)
    for i in range(len(drug)):
        score = SearchDrugGeneScore(drug.loc[i,'drug_claim_name'])
        tag = 'drug_claim_name'
        if score == []:
            score = SearchDrugGeneScore(drug.loc[i,'drug_claim_primary_name'])
            tag = 'drug_claim_primary_name'
            if score == []:
                score = SearchDrugGeneScore(drug.loc[i,'drug_name'])
                tag = 'drug_name'
                if score == []:
                    score = SearchDrugGeneScore(drug.loc[i,'drug_chembl_id'])
                    tag = 'drug_chembl_id'
                if score ==[]:
                    score = [['','','','','','']]
                    tag = ''
#         fi =['Drug_chembl_id|Drug_claim_primary_name','Drug_claim_primary_name','Drug_name','interactionTypes','geneName','sources','pmids','score','tag']
#         temp = fi
        for item in score:
            tem = [str(file.loc[i,'drug_claim_name']),str(file.loc[i,'drug_claim_primary_name']),
                   str(file.loc[i,'drug_name']),str(file.loc[i,'drug_chembl_id']),
                   item[1],item[2],item[3],item[4],item[5],tag]
            output = "\t".join(tem)
            output = output + '\n'
            f.writelines(output)

