#!/bin/bash
echo "Generating hardware info..."

# Get the hostname of the machine
DEVICE_NAME=$(hostname)
OUTPUT_FILE="${DEVICE_NAME}_hw.info"

#!/bin/bash
echo "Generating hardware info..."
{
    echo "CPU Info:"
    grep -m 1 'Model name' /proc/cpuinfo
    grep 'Physical CPU cores' /proc/cpuinfo | uniq
    echo "Logical processors: $(nproc --all)"
    
    echo ""
    echo "Cache Sizes:"
    echo -n "L1d Cache: "; cat /sys/devices/system/cpu/cpu0/cache/index0/size 2>/dev/null
    echo -n "L1i Cache: "; cat /sys/devices/system/cpu/cpu0/cache/index1/size 2>/dev/null
    echo -n "L2 Cache:  "; cat /sys/devices/system/cpu/cpu0/cache/index2/size 2>/dev/null
    echo -n "L3 Cache:  "; cat /sys/devices/system/cpu/cpu0/cache/index3/size 2>/dev/null

    echo ""
    echo "RAM Info:"
    grep MemTotal /proc/meminfo

    echo ""
    echo ""
    echo ""
    echo "====================="
    echo "All CPU info"
    echo "====================="
    cat /proc/cpuinfo

    echo ""
    echo "====================="
    echo "All memory info"
    echo "====================="
    cat /proc/meminfo
} > "$OUTPUT_FILE"

echo "Hardware info written to $OUTPUT_FILE"
