import pandas as pd
from collections import defaultdict
from rdkit import RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog('rdApp.*')

TRAIN_FRAC = 0.7   # 0.8 yapmak istersen buradan değiştir

df = pd.read_csv('data/egfr_clean.csv')

df['scaffold'] = df['smiles'].map(
    lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False))

groups = defaultdict(list)
for idx, scaf in zip(df.index, df['scaffold']):
    groups[scaf].append(idx)

# Büyük scaffold grupları önce train'e, sonra kalanlar test'e
ordered = sorted(groups.values(), key=lambda g: (-len(g), g[0]))
n_train_target = int(TRAIN_FRAC * len(df))

train_idx, test_idx = [], []
for g in ordered:
    if len(train_idx) + len(g) <= n_train_target:
        train_idx.extend(g)
    else:
        test_idx.extend(g)

train, test = df.loc[train_idx], df.loc[test_idx]

overlap = set(train['scaffold']) & set(test['scaffold'])
print(f'Toplam molekül: {len(df)}, benzersiz scaffold: {len(groups)}')
print(f'Train: {len(train)} ({len(train)/len(df):.0%}), Test: {len(test)} ({len(test)/len(df):.0%})')
print(f'Train/test ortak scaffold sayısı: {len(overlap)}  (0 olmalı)')
print(f"pChEMBL ortalaması -> train: {train['pchembl_median'].mean():.2f}, "
      f"test: {test['pchembl_median'].mean():.2f}")

train.to_csv('data/train.csv', index=False)
test.to_csv('data/test.csv', index=False)