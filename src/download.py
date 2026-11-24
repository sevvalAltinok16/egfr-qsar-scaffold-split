import pandas as pd
from chembl_webresource_client.new_client import new_client

activities = new_client.activity.filter(
    target_chembl_id='CHEMBL203',   # EGFR
    standard_type='IC50',
    pchembl_value__isnull=False,
).only(['molecule_chembl_id', 'canonical_smiles', 'standard_type',
        'standard_relation', 'standard_value', 'standard_units',
        'pchembl_value'])

df = pd.DataFrame(list(activities))
df.to_csv('data/egfr_raw.csv', index=False)
print(df.shape)
print(df.columns.tolist())