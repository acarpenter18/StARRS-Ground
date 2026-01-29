# Yan Li

# Libraries
from pymavlink import mavutil
import paramiko

# Set up
master = mavutil.mavlink_connection('/dev/ttyUSB0',baud=57600)
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('antpi', username='amc2918', password='')

# Main loop
while True:
    msg=master.recv_match()
    if msg:
        print(f"received message:{msg}")
        if 'start' in str(msg):
            ssh.exec_command('change_code_to_start_radar_py')
            print("radar starts")
        elif 'stop' in str(msg):
            ssh.exec_command('change_code_to_start_radar_py')
            print("radar stops")
        elif 'time' in str(msg):
            ssh.exec_command('uptime')
            print("system time:", stdout.read().decode())
