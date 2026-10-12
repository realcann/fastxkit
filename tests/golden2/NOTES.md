Davranış notları (referans araç gözlemleri):
- `sort` kopya kimlikli (aynı ID) girdide hata verir ve çıkış kodu 1 döner ("duplicated sequences found"); `sort -n/-l/-s` bu yüzden uniq.fa üzerinde.
- `grep -p desen` varsayılan olarak kayıt KİMLİĞİNE (ilk sözcük) bakar; `-n` tüm ada, `-s` diziye bakar.
- `-r` desenin düzenli ifade olduğunu, `-i` büyük/küçük harf duyarsızlığını belirtir.
