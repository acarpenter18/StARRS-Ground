#!/bin/bash
echo "Starting SAR Data Receiver"
echo "=========================="
echo ""
echo "Checking SSD mount..."
if mountpoint -q /mnt/ssd; then
    echo "SSD is mounted"
    df -h /mnt/ssd
else
    echo "SSD not mounted. Mounting now..."
    sudo mount /dev/sda1 /mnt/ssd
    if [ $? -eq 0 ]; then
        echo "SSD mounted successfully"
    else
        echo "Failed to mount SSD"
        echo "Using local storage instead"
        mkdir -p ~/sar_data
        ./rpi_1809 --path ~/sar_data --port 5555
        exit
    fi
fi

echo ""
echo "Starting receiver on port 5555..."
./rpi_1809 --path /mnt/ssd --port 5555
