import numpy as np
import pandas as pd
import pdb

def readin(fnName):
    df = pd.read_csv(fnName,sep='\t', index_col = 0, low_memory=False)
    df[df == '.'] = 0.0
    df[pd.isna(df)]   = 'None'
    df_numeric = df.convert_objects(convert_numeric=True)
    #df_numeric = pd.to_numeric(df)
    return df_numeric


