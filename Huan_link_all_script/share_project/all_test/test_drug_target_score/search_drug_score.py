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

os.chdir('/f/mulinlab/sinan/drug_reposition/DGIdb/')
file = pd.read_table('sort_unique_no_indication_drug_need_to_find_drug_target_score_left.txt',sep='\t')
file = file.loc[:,['Drug_chembl_id|Drug_claim_primary_name','Drug_claim_primary_name','Drug_name']]

with open('sort_unique_no_indication_drug_target_score_left.txt','w') as f:
    title = 'Drug_chembl_id|Drug_claim_primary_name\tDrug_claim_primary_name\tDrug_name\tinteractionTypes\tgeneName\tsources\tpmids\tscore\ttag\n'
    f.writelines(title)
    for i in range(len(file)):
        score = SearchDrugGeneScore(file.loc[i,'Drug_chembl_id|Drug_claim_primary_name'])
        tag = 'Drug_chembl_id|Drug_claim_primary_name'
        if score == []:
            score = SearchDrugGeneScore(file.loc[i,'Drug_claim_primary_name'])
            tag = 'Drug_claim_primary_name'
            if score == []:
                score = SearchDrugGeneScore(file.loc[i,'Drug_name'])
                tag = 'Drug_claim_primary_name'
                if score ==[]:
                    score = [['','','','','','']]
                    tag = ''
#         fi =['Drug_chembl_id|Drug_claim_primary_name','Drug_claim_primary_name','Drug_name','interactionTypes','geneName','sources','pmids','score','tag']
#         temp = fi
        for item in score:
            tem = [str(file.loc[i,'Drug_chembl_id|Drug_claim_primary_name']),str(file.loc[i,'Drug_claim_primary_name']),str(file.loc[i,'Drug_name']),
                   item[1],item[2],item[3],item[4],item[5],tag]
            output = "\t".join(tem)
            output = output + '\n'
            f.writelines(output)