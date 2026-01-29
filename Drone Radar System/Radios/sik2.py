# Yan Li
# Libraries
from pymavlink import mavutil
import paramiko

# ssh variables
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
sshparamiko = '/home/amc2918/sshparamiko.py'

# ssh details
hostname = 'antpi'
port = 22
username = 'amc2918'
password = '123456'

# Radio Set up
master = mavutil.mavlink_connection('/dev/ttyUSB0',baud=57600)

# ssh set up
client.connect(hostname, port, username, password)
ettus = f'python3 {sshparamiko}'
stdin, stdout, stderr = client.exec_command(ettus)

# Main loop
while True:
    msg=master.recv_match()
    if msg:
        print(f"received message:{msg}")
        if 'sud' in str(msg):
            print("sudo screen /dev/ttyUSB0 115200")
            print(stdout.read().decode())
