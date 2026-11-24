import os
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
import matplotlib.pyplot as plt

RDLogger.DisableLog('rdApp.*')
os.makedirs('reports', exist_ok=True)

RAW, OUT = 'data/egfr_raw.csv', 'data/egfr_clean.csv'
REPORT, FIG = 'reports/data_quality_report.txt', 'reports/pchembl_hist.png'

log = []
def step(name, df):
    log.append(f'{name}: {len(df)} kayıt')
    print(log[-1])

df = pd.read_csv(RAW, sep=None, engine='python')
step('Ham veri', df)

# --- 1. Aktivite filtreleri ---
df = df[df['Target ChEMBL ID'] == 'CHEMBL203']
step('Hedef = CHEMBL203', df)

df = df[df['Standard Type'] == 'IC50']
step('Tip = IC50', df)

rel = df['Standard Relation'].astype(str).str.strip("'\" ")
df = df[rel == '=']
step("Relation = '='", df)

df = df[df['Standard Units'] == 'nM']
step('Birim = nM', df)

df = df.copy()
df['pchembl'] = pd.to_numeric(df['pChEMBL Value'], errors='coerce')
df = df[df['pchembl'].notna()]
step('pChEMBL değeri var', df)

df = df[df['Data Validity Comment'].isna()]
step('Veri geçerlilik uyarısı yok', df)

dup = pd.to_numeric(df['Potential Duplicate'], errors='coerce').fillna(0)
df = df[dup == 0]
step('Potansiyel duplicate işaretli değil', df)

df = df[df['Assay Variant Mutation'].isna()]
step('Mutant olmayan (wild-type) assay', df)

df = df[df['Smiles'].notna()]
step('SMILES var', df)

# --- 2. Molekül standardizasyonu ---
n_salt = int(df['Smiles'].astype(str).str.contains(r'\.').sum())
chooser = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()

def standardize(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = chooser.choose(mol)
    mol = uncharger.uncharge(mol)
    return Chem.MolToSmiles(mol)

df['smiles_clean'] = df['Smiles'].astype(str).map(standardize)
n_invalid = int(df['smiles_clean'].isna().sum())
df = df[df['smiles_clean'].notna()]
step('Geçerli SMILES (standardizasyon sonrası)', df)

# --- 3. Duplicate'leri median ile birleştir ---
n_before = len(df)
g = df.groupby('smiles_clean')['pchembl']
out = g.agg(pchembl_median='median', n_measure='count', pchembl_std='std').reset_index()
out = out.rename(columns={'smiles_clean': 'smiles'})
n_noisy = int((out['pchembl_std'] > 1).sum())
step('Benzersiz molekül (median ile birleştirilmiş)', out)

out.to_csv(OUT, index=False)

# --- 4. Rapor ve grafik ---
lines = ['VERİ KALİTESİ RAPORU - EGFR (CHEMBL203)', '=' * 40] + log + [
    '',
    f'Tuz/solvent içeren kayıt (temizlenmeden önce): {n_salt}',
    f'Geçersiz SMILES nedeniyle elenen: {n_invalid}',
    f'Birleştirilen duplicate kayıt sayısı: {n_before - len(out)}',
    f'Ölçümleri arasında std > 1 olan molekül: {n_noisy}',
    f"pChEMBL medyan: {out['pchembl_median'].median():.2f}  "
    f"min: {out['pchembl_median'].min():.2f}  max: {out['pchembl_median'].max():.2f}",
]
with open(REPORT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n'.join(lines))

plt.hist(out['pchembl_median'], bins=30)
plt.xlabel('pChEMBL (IC50)'); plt.ylabel('Molekül sayısı')
plt.title('EGFR - temizlenmiş veri seti')
plt.savefig(FIG, dpi=150, bbox_inches='tight')