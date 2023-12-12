### What is the command to uncompress the file test.gz?

gzip -d test.gz `
gzip -r test.gz
gzip -c test.gz
gzip -v test.gz



### What is the command to set the following permission for a directory Biol3300?

User: read, write, execute

Group: read, execute

Others: read

chmod 745 Biol3300
chmod -R 754 Biol3300`
chmod 754 Biol3300
chmod -R 745 Biol3300

### The permission for a file are as follows "-rwxr-xr--". Who has execute permissions?

Users and groups `
Others
Groups
Users

### What option of gzip command can preserve the input files during compression or decompression?
-k
-p
-r
-v

### What command can you use to find the line number of the exon with the id "ENSMODE00000377263" from the gff3file?

grep -n ENSMODE00000377263 gff3file `     
grep -l ENSMODE00000377263 gff3file 
grep -c ENSMODE00000377263 gff3file 
grep -i ENSMODE00000377263 gff3file 


### What exon id is 2 lines before the exon with the id "ENSMODE00000377263" from the gff3file? 
ENSMODE00000370841 
ENSMODE00000315663 `
ENSMODE00000377263
ENSMODE00000321910


### How many exons are in the gff3file? 
25183 `
25184
25185
25186


### To create an archive named abc.tar consisting of 2 files, file1 and file2, which of the following command will be used?
tar -c abc.tar file1 file2
tar -cvf abc.tar file1 file2 `
tar -cv abc.tar file1 file2
tar -cvf file1 file2 abc.tar

### Which one of the following commands can add a new line "I like fruits" to the end of the fruit.txt file?
echo "I like fruits" > fruit.txt
echo "I like fruits" >> fruit.txt `
cat "I like fruits" > fruit.txt
cat "I like fruits" >> fruit.txt

### There is a gene sample file, gene.txt, as shown below:

aatctCGaaActaggAATaaggtactatg

What command could you use to convert the lower case a base to upper case A base?
tr 'A' 'a' gene.txt
tr 'a' 'A' < gene.txt `
tr 'a' 'A' gene.txt
tr 'A' 'a' > gene.txt

