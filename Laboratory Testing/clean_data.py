import pandas as pd

file_path = r'C:\Users\hp\OneDrive - Imperial College London\Desktop\output data\sweep_data.csv'


data = pd.read_csv(file_path)


expected_length = 501  # We have 501 points

# init
valid_sequences = []
current_sequence = []


for i in range(len(data)):
    if not current_sequence:
        current_sequence.append(data.iloc[i])
    elif data.iloc[i, 0] > data.iloc[i - 1, 0]:  # Still part of an increasing sequence
        current_sequence.append(data.iloc[i])
    else:  # Sequence ended
        if len(current_sequence) == expected_length:
            valid_sequences.extend(current_sequence)
        current_sequence = [data.iloc[i]]  # Start new sequence

if len(current_sequence) == expected_length:
    valid_sequences.extend(current_sequence)

processed_data = pd.DataFrame(valid_sequences, columns=data.columns)

processed_data['Frequency'] = processed_data['Frequency'].astype('Int64')

processed_data.to_csv(r'C:\Users\hp\OneDrive - Imperial College London\Desktop\output data\corrected_sweep_data.csv', index=False)

