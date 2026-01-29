#!/bin/bash
# optimise_streaming.sh - Streaming optimisation for E312 FMCW data collection

echo "======================================"
echo "   E312 Optimisation Script  "
echo "======================================"
echo ""

# Function to check if command succeeded
check_status() {
    if [ $? -eq 0 ]; then
        echo "  ✓ $1"
    else
        echo "  ✗ $1 (failed, but continuing...)"
    fi
}

# 1. CPU Optimisation
echo "1. CPU Optimisation:"

# Disable CPU idle states for lower latency
if [ -d "/sys/devices/system/cpu/cpu0/cpuidle" ]; then
    for state in /sys/devices/system/cpu/cpu*/cpuidle/state*/disable; do
        echo 1 > $state 2>/dev/null
    done
    check_status "CPU idle states disabled"
fi

# Check current CPU idle governor (read-only, just for info)
if [ -f "/sys/devices/system/cpu/cpuidle/current_governor_ro" ]; then
    gov=$(cat /sys/devices/system/cpu/cpuidle/current_governor_ro)
    echo "  ℹ CPU idle governor: $gov (read-only)"
fi

# Set CPU affinity for interrupts to CPU 0 (leave CPU 1 for application)
echo 1 > /proc/irq/default_smp_affinity 2>/dev/null
check_status "IRQ affinity set to CPU 0"

# 2. Memory Optimisation
echo ""
echo "2. Memory Optimisation:"

# Disable swap
echo 0 > /proc/sys/vm/swappiness 2>/dev/null
check_status "Swappiness set to 0"

# Clear caches
sync
echo 3 > /proc/sys/vm/drop_caches 2>/dev/null
check_status "Page cache cleared"

# Optimize memory allocation
echo 1 > /proc/sys/vm/overcommit_memory 2>/dev/null
check_status "Memory overcommit enabled"

# 3. Network Buffer Optimisation
echo ""
echo "3. Network Optimisation:"

# Increase network buffers
sysctl -w net.core.rmem_max=134217728 >/dev/null 2>&1
check_status "Receive buffer increased to 128MB"

sysctl -w net.core.wmem_max=134217728 >/dev/null 2>&1
check_status "Send buffer increased to 128MB"

sysctl -w net.core.netdev_max_backlog=5000 >/dev/null 2>&1
check_status "Network backlog increased"

# TCP optimizations for ZMQ
sysctl -w net.ipv4.tcp_nodelay=1 >/dev/null 2>&1
check_status "TCP nodelay enabled"

# Find active network interface and optimize MTU
echo ""
echo "  Network interfaces:"
for iface in $(ls /sys/class/net/); do
    if [ "$iface" != "lo" ]; then
        current_mtu=$(cat /sys/class/net/$iface/mtu 2>/dev/null)
        if [ -n "$current_mtu" ]; then
            # Try MTU 4000 (more realistic than 9000)
            ip link set $iface mtu 4000 2>/dev/null
            if [ $? -eq 0 ]; then
                echo "    $iface: MTU $current_mtu → 4000"
            else
                # Try 2000 as fallback
                ip link set $iface mtu 2000 2>/dev/null
                if [ $? -eq 0 ]; then
                    echo "    $iface: MTU $current_mtu → 2000"
                else
                    echo "    $iface: MTU $current_mtu (unchanged)"
                fi
            fi
        fi
    fi
done

# 4. USB Optimisation
echo ""
echo "4. USB Optimisation:"

if [ -f "/sys/module/usbcore/parameters/autosuspend" ]; then
    echo -1 > /sys/module/usbcore/parameters/autosuspend 2>/dev/null
    check_status "USB autosuspend disabled"
fi

# 5. Process Scheduling
echo ""
echo "5. Process Scheduling:"

# Enable real-time scheduling
if [ -f "/proc/sys/kernel/sched_rt_runtime_us" ]; then
    echo -1 > /proc/sys/kernel/sched_rt_runtime_us 2>/dev/null
    check_status "Real-time scheduling unlimited"
fi

# Reduce scheduler latency
echo 1000000 > /proc/sys/kernel/sched_min_granularity_ns 2>/dev/null
check_status "Scheduler granularity optimised"

# 6. Stop Unnecessary Services
echo ""
echo "6. Stopping unnecessary services:"

services="nginx apache2 bluetooth avahi-daemon cups ntp chronyd"
for service in $services; do
    if systemctl is-active --quiet $service 2>/dev/null; then
        systemctl stop $service 2>/dev/null
        check_status "Stopped $service"
    fi
done

# 7. Kernel Module Optimisation
echo ""
echo "7. Kernel Optimisation:"

# Disable kernel debug features if present
echo 0 > /proc/sys/kernel/ftrace_enabled 2>/dev/null
check_status "Ftrace disabled"

echo 0 > /proc/sys/kernel/nmi_watchdog 2>/dev/null
check_status "NMI watchdog disabled"

# 8. File System Optimisation
echo ""
echo "8. File System:"

# Use noatime for reduced disk writes
mount -o remount,noatime / 2>/dev/null
check_status "Filesystem mounted with noatime"

# 9. Power Management
echo ""
echo "9. Power Management:"

# Disable frequency scaling
for gov in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
    if [ -f "$gov" ]; then
        echo performance > $gov 2>/dev/null
        check_status "CPU frequency governor set to performance"
        break
    fi
done

# Zynq-specific power management
if [ -d "/sys/devices/soc0" ]; then
    echo "  ℹ Zynq SoC detected"
fi

# 10. System Status Report
echo ""
echo "======================================"
echo "System Status:"
echo "======================================"

# CPU info
echo ""
echo "CPU Information:"
nproc=$(nproc)
echo "  Cores: $nproc"
grep "BogoMIPS\|MHz" /proc/cpuinfo | head -2

# Memory
echo ""
echo "Memory Status:"
free -h | grep -E "^Mem|^Swap" | while read line; do echo "  $line"; done

# Temperature
echo ""
echo "Temperature:"
for zone in /sys/class/thermal/thermal_zone*/temp; do
    if [ -f "$zone" ]; then
        temp=$(cat $zone)
        temp_c=$((temp / 1000))
        echo "  $(basename $(dirname $zone)): ${temp_c}°C"
    fi
done

# Network
echo ""
echo "Network Configuration:"
ip -br link show | grep -v "lo" | while read line; do echo "  $line"; done

# Process priority hint
echo ""
echo "======================================"
echo "Optimisation Complete!"
echo "======================================"

# Create a flag file to indicate optimisation was run
touch /tmp/e312_optimized
echo "Optimisation flag set: /tmp/e312_optimised"