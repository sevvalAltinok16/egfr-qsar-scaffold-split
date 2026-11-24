from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

for ad, smi in [('etanol', 'CCO'), ('propanol', 'CCCO'), ('benzen', 'c1ccccc1')]:
    fp = gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(smi))
    print(ad, '| uzunluk:', len(fp), '| 1 olan kutu sayısı:', fp.sum(),
          '| kutular:', fp.nonzero()[0])