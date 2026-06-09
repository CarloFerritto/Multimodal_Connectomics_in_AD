import pandas as pd
import numpy as np
from neuroHarmonize import harmonizationLearn, harmonizationApply


# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================
print("📂 Loading and Preparing Data...")
data_dir = "./data/"
meta_df = pd.read_excel(data_dir + 'dataset_summary_ALL.xlsx',sheet_name='ALL')
df = pd.read_csv(data_dir + 'func.csv')
df_other = pd.read_csv(data_dir + 'func.csv')
ref_df = pd.read_excel(data_dir + 'subjects_harmonization.xlsx')
df_transposed = df.T
n_features = df_transposed.shape[1]
df_transposed.columns = [f"Feature_{i+1}" for i in range(n_features)]

df_other_transposed = df_other.T
n_features_other = df_other_transposed.shape[1]
df_other_transposed.columns = [f"Feature_{i+1}" for i in range(n_features_other)]
df_other_transposed.insert(0, "Subject", df_other_transposed.index)
df_other_transposed = df_other_transposed.reset_index(drop=True)
# Add the "Subject" column with original indices
df_transposed.insert(0, "Subject", df_transposed.index)

# Reset index
df_transposed = df_transposed.reset_index(drop=True)

print(df_transposed.head())


# Merge & Filter Common Subjects
common_ids = set(meta_df['Subject']).intersection(df_transposed['Subject'])
meta = meta_df[meta_df['Subject'].isin(common_ids)].sort_values('Subject').reset_index(drop=True)
feat = df_transposed[df_transposed['Subject'].isin(common_ids)].sort_values('Subject').reset_index(drop=True)
feat_other = df_other_transposed[df_other_transposed['Subject'].isin(common_ids)].sort_values('Subject').reset_index(drop=True)
assert all(meta['Subject'] == feat['Subject'])
assert all(meta['Subject'] == feat_other['Subject'])

# Diagnosis Handling
if 'Diagnosis_amyloid' not in meta.columns:
    meta['Diagnosis_amyloid'] = meta['Diagnosis']

# Batch Filter (< 2 CN-)
batch_counts = meta[meta['Diagnosis_amyloid'] == 'CN-'].groupby('Batch_variable_simplified').size()
valid_sites = batch_counts[batch_counts >= 2].index.tolist()
meta_clean = meta[meta['Batch_variable_simplified'].isin(valid_sites)].reset_index(drop=True)
feat_clean = feat[feat['Subject'].isin(meta_clean['Subject'])].reset_index(drop=True)
feat_other_clean = feat_other[feat_other['Subject'].isin(meta_clean['Subject'])].reset_index(drop=True)
print(f"✅ Subjects ready for analysis: {len(meta_clean)}")

# Data Matrices
data_matrix = feat_clean.drop(columns=['Subject']).values
data_other_matrix = feat_other_clean.drop(columns=['Subject']).values  # not harmonized
features_names = feat_clean.drop(columns=['Subject']).columns.tolist()
features_other_names = feat_other_clean.drop(columns=['Subject']).columns.tolist()
covars = meta_clean[['Batch_variable_simplified', 'Age', 'Sex', 'Diagnosis_amyloid']].copy()
# ==============================================================================
# 2. SAVING SCENARIOS
# ==============================================================================

# --- SCENARIO 1: Reference-Based ---
print("\n🚀 Running SCENARIO 1 (Reference-Based)...")
is_ref = meta_clean['Subject'].isin(ref_df['Subject'])
idx_train = meta_clean.index[is_ref].tolist()

if len(idx_train) == 0:
    raise ValueError("No Reference subjects found!")

data_train = data_matrix[idx_train, :]
site_series = 'Site_' + meta_clean['Batch_variable_simplified'].astype(str)
sex_series = meta_clean['Sex'].astype('category').cat.codes
# Base covars DataFrame (for Scenario 1 and 2)
# neuroHarmonize expects the column "SITE" in the covars dataframe
covars_base = pd.DataFrame({
    'SITE': site_series.values,
    'AGE': meta_clean['Age'].values,
    'SEX': sex_series.values

})
covars_train = covars_base.iloc[idx_train].reset_index(drop=True)

# Learn (on Reference only)
# smooth_terms=[] indicates linear model for age (standard ComBat)
model_s1, _ = harmonizationLearn(data_train, covars_train, smooth_terms=[])

# Apply (on ALL subjects)
# Important: covars_base must have the same columns as covars_train
data_s1_adj = harmonizationApply(data_matrix, covars_base, model_s1)
data_s1_adj_df = pd.DataFrame(data_s1_adj, columns=features_names)
data_s1_adj_df.insert(0, "Subject", meta_clean['Subject'].values)

# Transpose back to original format
data_s1 = data_s1_adj_df.set_index("Subject").T

# Save S1
data_s1.to_csv(data_dir + "func_combat.csv", header=True, index=False)

data_other_df = pd.DataFrame(data_other_matrix, columns=features_other_names)
data_other_df.insert(0, "Subject", meta_clean['Subject'].values)
data_other = data_other_df.set_index("Subject").T
data_other.to_csv(data_dir + "func_new.csv", header=True, index=False)