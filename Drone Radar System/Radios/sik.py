# Yan Li

from pymavlink import mavutil
master=mavutil.mavlink_connection('/dev/ttyUSB0',baud=57600)
while True:
    msg=master.recv_match()
    if msg:
        print(f"received message:{msg}")
        session_name_or_id = "1319"
        text_to_send = "sudo xxx"
        text_to_send_enter = "./New 5.4e9 5e6 50 20\n"
        """
        print(f"screen -S {session_name_or_id} -X stuff \"{text_to_send}\"")
        import os
        command = f"screen -S {session_name_or_id} -X stuff \"{text_to_send}\""
        os.system(command)
        """
        import subprocess
        command = ["screen", "-S", session_name_or_id, "-X", "stuff", text_to_send_enter]
        subprocess.run(command, check=True)
        
