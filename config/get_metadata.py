import pandas as pd
import re
from pathlib import Path


def basecall_version(x):
    r = re.search('(\d+\.\d+\.\d+)', str(x))
    if r:
        return r.group(1)


def basecall_model(x):
    r = re.search('\.(hac|fast)\.', str(x), flags=re.I)
    if r:
        return r.group(1)


root = Path(__file__).parents[1]
df = pd.read_csv(root/'datasets.csv', comment='#')
print(df)

# FAO47403:
#     sample_name:
#     genome:
#     basecall:
#     model:

for idx, dfx in df.iterrows():
    print(
        f'''
        {dfx.run_id}:
            method:
            genome: {dfx.dataset_id}
            basecall_version: {basecall_version(dfx.fastq_path)}
            basecall_model: {basecall_model(dfx.fastq_path)}
        '''
    )
    # print(dfx.fastq_path.split('/')[-1])
    # print(basecall_version(dfx.fastq_path))
    # print(basecall_model(dfx.fastq_path))
