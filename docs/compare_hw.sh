#!/bin/bash

# Usage: ./compare_hardware.sh sysinfo1.txt sysinfo2.txt

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 file1 file2"
  exit 1
fi

FILE1="$1"
FILE2="$2"

# Check both files exist
for f in "$FILE1" "$FILE2"; do
  if [ ! -f "$f" ]; then
    echo "Error: File '$f' not found."
    exit 1
  fi
done

# Function to extract a value from a file
get_value() {
  grep -m 1 "$1" "$2" | awk -F ':' '{gsub(/^[ \t]+/, "", $2); print $2}'
}

# Features to compare
declare -A LABELS=(
  ["model name"]="CPU Model"
  ["cpu cores"]="Physical cores"
  ["Logical processors"]="Logical processors"
  ["L1d Cache"]="L1d Cache"
  ["L1i Cache"]="L1i Cache"
  ["L2 Cache"]="L2 Cache"
  ["L3 Cache"]="L3 Cache"
  ["MemTotal"]="Total RAM"
)

echo "=== Hardware Comparison ==="
printf "%-20s | %-45s | %-45s\n" "Feature" "$(basename "$FILE1" _hw.info)" "$(basename "$FILE2" _hw.info)"
printf -- "---------------------|-----------------------------------------------|--------------------------------------------------\n"

ordered_keys=("model name" "cpu cores" "Logical processors" "L1d Cache" "L1i Cache" "L2 Cache" "L3 Cache" "MemTotal")
for key in "${ordered_keys[@]}"; do
    val1=$(get_value "$key" "$FILE1")
    val2=$(get_value "$key" "$FILE2")
    printf "%-20s | %-45s | %-45s\n" "${LABELS[$key]}" "$val1" "$val2"
done
