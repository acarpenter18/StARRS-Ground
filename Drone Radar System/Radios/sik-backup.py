# Yan Li

from pymavlink import mavutil
master=mavutil.mavlink_connection('/dev/ttyUSB0',baud=57600)
while True:
    msg=master.recv_match()
    if msg:
        print(f"received message:{msg}")
