#!/usr/bin/env bash
# Referans çıktıları üretir (yalnızca sentetik tests/data üzerinde). Kullanım: make_golden.sh <referans-ikili>
# Çıktılar golden/*.out olarak kaydedilir; testler ikiliye ihtiyaç duymaz.
set -u
SK="$1"; cd "$(dirname "$0")/../data"; OUT=../golden
run() { name="$1"; shift; "$SK" "$@" 2>/dev/null > "$OUT/$name.out"; }
run stats_basic        stats small.fa small.fq small.fq.gz
run stats_all_fq       stats -a small.fq
run stats_all_fa       stats -a small.fa
run seq_rc_fa          seq -r -p small.fa
run seq_upper_fa       seq -u small.fa
run seq_lower_fa       seq -l small.fa
run seq_name_fa        seq -n small.fa
run seq_id_fa          seq -n -i small.fa
run seq_seq_fa         seq -s small.fa
run seq_w0_fa          seq -w 0 small.fa
run seq_w30_fa         seq -w 30 small.fa
run seq_minlen_fa      seq -m 60 small.fa
run seq_maxlen_fq      seq -M 100 small.fq
run seq_rc_fq          seq -r -p small.fq
run head_fa            head -n 2 small.fa
run head_fq            head -n 3 small.fq
run range_fa           range -r 2:4 small.fa
run range_neg_fa       range -r -2:-1 small.fa
run range_fq           range -r 1:2 small.fq
run fx2tab_fa          fx2tab small.fa
run fx2tab_fq          fx2tab small.fq
run fx2tab_lg_fq       fx2tab -l -g small.fq
run fx2tab_name_fq     fx2tab -n small.fq
run fx2tab_id_fa       fx2tab -n -i small.fa
run fx2tab_q_fq        fx2tab -q small.fq
run fx2tab_lg_fa       fx2tab -l -g -B AT small.fa
"$SK" fx2tab small.fa 2>/dev/null | "$SK" tab2fx 2>/dev/null > "$OUT/tab2fx_fa.out"
"$SK" fx2tab small.fq 2>/dev/null | "$SK" tab2fx 2>/dev/null > "$OUT/tab2fx_fq.out"
