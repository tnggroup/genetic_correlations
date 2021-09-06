import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Remove duplicates and add N columns')
parser.add_argument('-i','--input', help='Input sumstats (will be overriden)', required=True)
parser.add_argument('-a','--ntotal', help='Total sample size')
parser.add_argument('-b','--ncas', help='Sample size cases')
parser.add_argument('-c','--ncon', help='Sample size controls')

args = vars(parser.parse_args())

df = pd.read_csv(args['input'],sep="\t")

df.dropna(inplace=True)

checkP =  'P' in df
checkFreq = 'FREQ' in df
checkFreqCases = 'FREQCASES' in df
checkInfo = 'INFO' in df
checkOr = 'OR' in df

def findfield(df):
    truefalse=True
    if checkP:
        truefalse = (df['P'] >= 0.0) and (df['P'] <= 1.0)
    else:
        truefalse = False
    if checkFreq: 
        truefalse = truefalse and (df['FREQ'] >= 0.005) and (df['FREQ'] <= 0.995)
    if checkFreqCases: 
        truefalse = truefalse and (df['FREQCASES'] >= 0.005) and (df['FREQCASES'] <= 0.995)
    if checkInfo: 
        truefalse = truefalse and (df['INFO'] >= 0.6)
    if checkOr: 
        truefalse = truefalse and (df['OR'] <= 10000)
        return truefalse

df['truefalse'] = df.apply(lambda x: findfield(x), axis=1)

df = df.loc[df['truefalse']]

df.drop(columns=['truefalse'],inplace=True) 

df.drop_duplicates(['SNP'], keep=False, inplace=True)

if ((args['ntotal'] is not None) and (args['ntotal']!='')) and ((args['ncas'] is not None) and (args['ncas']!='')) and ((args['ncon'] is not None) and (args['ncon']!='')) and ('N' not in df) and ('Ncas' not in df) and ('Ncon' not in df):
    df['N'] = args['ntotal']
    df['Ncas'] = args['ncas']
    df['Ncon'] = args['ncon']
elif (args['ntotal'] is not None and (args['ntotal']!='')) and ('N' not in df):
    df['N'] = args['ntotal']
 
df.to_csv(args['input'], compression='gzip',sep="\t",index=False)

