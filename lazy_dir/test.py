import pandas as pd

df = pd.read_csv('preprocessor_rejects.csv')

print(df['Tool'].value_counts())