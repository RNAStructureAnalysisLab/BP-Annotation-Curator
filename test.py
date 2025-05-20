import os
import pandas as pd

path = os.path.join('Data', 'Raw', 'AnnotationTools', 'Annotations', '4V9F_CL.csv')
df = pd.read_csv(path, dtype={'residue1': str, 'residue2': str})
print(df['residue1'].dtype)

print(df.head(50))