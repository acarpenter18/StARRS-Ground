import time
import math
import pandas as pd
import os
from libreVNA import libreVNA
from tqdm import tqdm

# Create control interface
vna = libreVNA('localhost', 19542)  # SCPI server port 19542

# Quick connection check
print(vna.query("*IDN?"))

# Ensure connection to the device
dev = vna.query(":DEV:CONN?")
if dev == "Not connected":
    print("Not connected to any device, aborting")
    exit(-1)
else:
    print("Connected to " + dev)

# Stop current data acquisition
vna.cmd(":VNA:ACQ:STOP")

# Switch to VNA mode and set scan parameters
print("Setting up the sweep...")
vna.cmd(":DEV:MODE VNA")
vna.cmd(":VNA:SWEEP FREQUENCY")
vna.cmd(":VNA:STIM:LVL -10")
vna.cmd(":VNA:ACQ:IFBW 100")
vna.cmd(":VNA:ACQ:AVG 1")
vna.cmd(":VNA:ACQ:POINTS 501")
vna.cmd(":VNA:FREQ:START 1000000")
vna.cmd(":VNA:FREQ:STOP 6000000000")  # Modify the frequency stop point to 6000 MHz


# Define callback function to handle each data point
all_sweep_data = []
sweep_data = []
sweepComplete = False

def callback(data):
    global sweepComplete, sweep_data
    # print(f"Callback received data point: {data}")
    sweep_data.append(data)
    if data["pointNum"] == 500:
        # If this is the last point
        # print("Last data point received, stopping callback.")
        sweepComplete = True

num_sweeps = 5  # Set the number of sweeps to 5
for i in tqdm(range(num_sweeps)):
    sweepComplete = False
    sweep_data = []  # Clear current sweep data
    # print(f"Stopping previous sweep if any...")
    vna.cmd(":VNA:ACQ:STOP")  # Ensure the previous sweep is stopped
    vna.remove_live_callback(19000, callback)  # Ensure previous callback is removed
    # print(f"Starting sweep {i+1}...")
    vna.cmd(":VNA:ACQ:RUN")
    vna.add_live_callback(19000, callback)  # Real-time data stream port 19000
    

    while not sweepComplete:
        time.sleep(0.1)

    all_sweep_data.append(sweep_data)  # Save current sweep data
    # print(f"Sweep {i+1} complete, total data points collected: {len(sweep_data)}")
vna.cmd(":VNA:ACQ:STOP")
time.sleep(0.1)
print(f"Total sweeps completed: {num_sweeps}")

# Process or save all sweep data
all_processed_data = []
for sweep_num, sweep in enumerate(all_sweep_data):
    # print(f"Sweep {sweep_num+1} data:")
    for data_point in sweep:
        # Print data points to check their structure
        # print(f"Data point: {data_point}")
        
        frequency = data_point["frequency"]
        measurements = data_point["measurements"]
        
        s11_real = measurements["S11"].real if "S11" in measurements else 0
        s11_imag = measurements["S11"].imag if "S11" in measurements else 0
        
        s12_real = measurements["S12"].real if "S12" in measurements else 0
        s12_imag = measurements["S12"].imag if "S12" in measurements else 0
        
        s21_real = measurements["S21"].real if "S21" in measurements else 0
        s21_imag = measurements["S21"].imag if "S21" in measurements else 0
        
        s22_real = measurements["S22"].real if "S22" in measurements else 0
        s22_imag = measurements["S22"].imag if "S22" in measurements else 0

        all_processed_data.append([frequency, s11_real, s11_imag, s12_real, s12_imag, s21_real, s21_imag, s22_real, s22_imag])

print("All sweeps completed")

# Print current working directory
current_directory = os.getcwd()
print(f"Current working directory: {current_directory}")

# Determine file save path and avoid filename conflicts
output_file_path = os.path.join(current_directory, "sweep_data.xlsx")
counter = 1
while os.path.exists(output_file_path):
    output_file_path = os.path.join(current_directory, f"sweep_data_{counter}.xlsx")
    counter += 1

# Convert processed data to DataFrame and save to Excel file
df = pd.DataFrame(all_processed_data, columns=["Frequency", "S11_Real", "S11_Imag", "S12_Real", "S12_Imag", "S21_Real", "S21_Imag", "S22_Real", "S22_Imag"])
df.to_excel(output_file_path, index=False)
print(f"Data has been saved to {output_file_path}")
