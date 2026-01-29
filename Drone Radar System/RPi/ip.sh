#Set RPi IP to match Ettus for data streaming

ip addr flush dev eth0
ip addr add 192.168.1.100/24 dev eth0
ip link set eth0 up
ip addr
