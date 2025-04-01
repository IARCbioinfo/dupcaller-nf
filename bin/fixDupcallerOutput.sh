#!/bin/bash

for vcf in "$@"
do

out=${vcf/_calls.vcf/_fixed_calls.vcf}

first_format_num=$(grep -n -m 1 '##FORMAT' "$vcf" | cut -d : -f 1)
sed "${first_format_num}a##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" "$vcf" > "$out"
#sed -ri 's|(AC:)|GT:\1|' "$out"
sed -E 's/(AC:RC:DP[[:space:]]+)([^:]+:[^:]+:[^:]+)[[:space:]]+([^:]+:[^:]+:[^:]+)/GT:\10\/1:\2\t0\/0:\3/' "$vcf" > "$out"

gzip -f $out

done
