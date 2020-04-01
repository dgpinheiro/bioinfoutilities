#!/bin/sh

# Sugerido em http://www.ncbi.nlm.nih.gov/books/NBK47527/
# -l (maximum bandwidth of request, try 200M and go up from there) - Megabits/second
# -i <private key file> (http://pt.wikipedia.org/wiki/Criptografia_de_chave_p%C3%BAblica)
# -r recursive copy
# -Q (for adaptive flow control) - needed for disk throttling! (ajustar o fluxo de dados ao gravar no disco)
# -T to disable encryption
# -k1 enable resume of failed transfers
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -QTr -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP021/SRP021491/ .
