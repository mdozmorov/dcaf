#!/bin/bash

cut -b 1-3 $1 --complement | sort -k1 -k2,3 -n | python count_coverage.py
