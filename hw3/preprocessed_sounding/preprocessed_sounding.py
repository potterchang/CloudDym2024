import pandas as pd

# Load the CSV file
file_path = "no_2305_20210402_0559_L4.csv"
df = pd.read_csv(file_path)

# Extract relevant columns
df_processed = df[['P', 'T_new', 'q_new', 'U', 'V']].copy()

# Convert temperature from Celsius to Kelvin
df_processed['T_new'] = df_processed['T_new'] + 273.15
df_processed['U'] = 0
df_processed['V'] = 0

# Rename the columns to match the required format
df_processed.columns = ['P[hPa]', 'T[K]', 'Qv[g/kg]', 'U[m/s]', 'V[m/s]']

# Save the processed data to a text file
output_file_path = 'profile_sounding_QC.txt'
df_processed.to_csv(output_file_path, sep=' ', index=False, float_format='%.4f')

