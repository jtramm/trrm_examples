#!/bin/bash

awk 'BEGIN {
  max_percent = 0
} {
  if ($3 == "+/-") {
    value = $2
    std_dev = $4
    if (value + 0 != 0) {
      percent = (std_dev / value) * 100
      if (percent > max_percent) {
        max_percent = percent
      }
      printf("%s (%.6f%%)\n", $0, percent)
    } else {
      printf("%s (undefined%%)\n", $0)
    }
  } else {
    print $0
  }
} END {
  printf("\nMaximum percentage standard deviation: %.6f%%\n", max_percent)
}'

