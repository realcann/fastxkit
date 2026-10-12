# fastxkit: plan (2026-10-12)

Amaç: seqkit tarzı, yalnızca standart kütüphane kullanan Python FASTA/FASTQ aracı; değer ekseni: ham FASTQ'dan bağımsız QC sayıları ve olgu.json uyumlu çıktı.
Kaynak rapor: yazarın kişisel fikir notları (depoda yer almaz).

## Kararlar (Ruling)
- Ad `fastxkit` (marka karışıklığı olmasın). Orijinal kod okunmaz/kopyalanmaz; seqkit ikilisi yalnızca sentetik küçük dosyalarda doğruluk referansı (ağsız komutlar).
- Kapsam 21 komut, 3 dalga. Dışarıda: zstd, lz4, BAM, watch, scat, grafik.
- Hız hedefi yok.
- Her dalga: temiz-kod-yazici -> doğrulama (pytest + altın karşılaştırma) -> düzeltme.

## Dalgalar
1. G/Ç (gz/bz2/xz/stdin), ayrıştırıcı; stats, seq, head, range, sample, shuffle, fx2tab, tab2fx
2. grep, locate, subseq, faidx, sort, rmdup, rename, replace
3. fq2fa, split, sum, qc (olgu uyumlu JSON; yapıldı), crosscheck (yapıldı)

## Başarı ölçütü
seqkit ile birebir eşleşen çıktılar (altın testler) + qc'nin olgu.json demosu.
