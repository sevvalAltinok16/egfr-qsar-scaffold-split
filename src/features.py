import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog('rdApp.*')
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def featurize(smiles_list):
    return np.array([gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))
                     for s in smiles_list], dtype=np.uint8)

for name in ['train', 'test']:
    df = pd.read_csv(f'data/{name}.csv')
    X = featurize(df['smiles'])        # girdi: 0/1 vektörleri
    y = df['pchembl_median'].values    # hedef: tahmin edilecek aktivite
    np.save(f'data/X_{name}.npy', X)
    np.save(f'data/y_{name}.npy', y)
    print(name, 'X:', X.shape, 'y:', y.shape)