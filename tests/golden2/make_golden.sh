#!/usr/bin/env bash
# Dalga 2/3 referans çıktıları (yalnızca sentetik tests/data). Kullanım: make_golden.sh <referans-ikili>
set -u
SK="$1"; cd "$(dirname "$0")/../data"; OUT=../golden2
run() { name="$1"; shift; timeout 20 "$SK" "$@" 2>/dev/null </dev/null > "$OUT/$name.out"; }
run grep_id_fa         grep -p "seq2" small.fa
run grep_name_re_fa    grep -n -r -p "seq[24]" small.fa
run grep_invert_fa     grep -n -v -r -p "seq[24]" small.fa
run grep_seq_fa        grep -s -i -p "ACGT" dup.fa
run grep_seq_regex_fa  grep -s -r -p "A{3}G{3}" dup.fa
run grep_desc_fa       grep -n -r -p "record [13]" small.fa
run grep_seq_fq        grep -s -p "GGGCCC" motif.fq
run locate_fa          locate -p "ACGT" dup.fa
run locate_nonover_fa  locate -P -p "ACGT" dup.fa
run locate_regex_fq    locate -r -p "A+G+" motif.fq
run subseq_region_fa   subseq -r 2:10 small.fa
run subseq_neg_fa      subseq -r -5:-1 small.fa
run subseq_flank_fa    subseq -r 11:20 small.fa
run sort_name_fa       sort -n uniq.fa
run sort_len_fa        sort -l uniq.fa
run sort_len_rev_fa    sort -l -r uniq.fa
run sort_seq_fa        sort -s uniq.fa
run rmdup_seq_fa       rmdup -s dup.fa
run rmdup_seq_i_fa     rmdup -s -i dup.fa
run rmdup_name_fa      rmdup -n dup.fa
run rmdup_id_fa        rmdup dup.fa
run rename_fa          rename dup.fa
run replace_name_fa    replace -p "seq" -r "S" small.fa
run replace_re_name_fa replace -p "^seq(\\d)" -r "id_\$1" small.fa
run replace_seq_fa     replace -s -p "N" -r "-" small.fa
run replace_pat_fq     replace -p "^r" -r "read_" motif.fq
cp small.fa /tmp/_sf.fa; timeout 20 "$SK" faidx /tmp/_sf.fa 2>/dev/null </dev/null; cp /tmp/_sf.fa.fai "$OUT/faidx_small_fa.fai"
timeout 20 "$SK" faidx /tmp/_sf.fa seq1 seq5:10-20 2>/dev/null </dev/null > "$OUT/faidx_extract.out"
