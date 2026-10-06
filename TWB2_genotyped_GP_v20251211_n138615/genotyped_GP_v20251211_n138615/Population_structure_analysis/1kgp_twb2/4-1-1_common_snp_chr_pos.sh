#!/bin/bash

set -euo pipefail

mkdir -p chr_pos

cut -f2 ../chr_pos/1kgp_pos.bim | sort -u > 1kgp.pos
cut -f2 ../chr_pos/twb2_pos.bim  | sort -u > twb2.pos

comm -12 1kgp.pos twb2.pos > 1k_twb2_common.pos

# rm tmp.pos

echo "Common CHR:POS:"
wc -l 1k_twb2_common.pos