### If you use the command "ls *at" what files will be listed?
cat.txt, bat.txt, mat.txt, beat.txt
cat, bat, mat, sat
boat, beat, seat, meat
Both B and C

### If you use the command “ls ?at” what files will be listed?
cat.txt, bat.txt, mat.txt, beat.txt
cat, bat, mat, sat
boat, beat, seat, meat
Both B and C

### Here is an employee with salary and working years table, employee.txt

##### Manager 5000 3

##### Clerk 4000 5

##### Peon 4500 2

##### Director 9000 10

##### Guard 3000 2

### Which command can sort the table by their salary in numerical descending order? 
sort -k2 -n employee.txt
sort -n2 -r employee.txt 
sort -k2 employee.txt 
sort -k2 -n -r employee.txt 


### Here is an employee with salary and working years table, employee.txt

##### Manager 5000 3

##### Clerk 4000 5

##### Peon 4500 2

##### Director 9000 10

##### Guard 3000 4

### Which command can sort the table firstly by their working years and then by their salary in both numerical descending order? 
sort -k3nr -k2n employee.txt 
sort -k3 -nr employee.txt
sort -k3nr -k2nr employee.txt
sort -k2nr -k3nr employee.txt


### Here is an employee with salary and working years table, employee.txt

##### Manager 5000 3

##### Clerk 4000 5

##### Peon 4500 2

##### Director 9000 10

##### Guard 3000 4

### Which command can print out a new file with employee and their working years?
cut -f1 -f3 employee.txt 
cut -f1,3 employee.txt
cut -f1 -f3 -d"," employee.txt
cut -f1,3 -d" " employee.txt

### Which option is used with sort command for removing repeated lines?
-u
-t
-n
-a

### What is the default delimiter used by the cut command for cutting fields?
underscore
space
tab
comma


### Whatever we have cut using the cut command can be pasted back to the source file using the paste command but vertically.
True `
False

### In order to extract the first 20 lines of Monodelphis_domestica.ASM229v1.99.primary_assembly.6.gff3, and what command can you use?
head Monodelphis_domestica.ASM229v1.99.primary_assembly.6.gff3
tail Monodelphis_domestica.ASM229v1.99.primary_assembly.6.gff3
head 20 Monodelphis_domestica.ASM229v1.99.primary_assembly.6.gff3
head -n20 Monodelphis_domestica.ASM229v1.99.primary_assembly.6.gff3

### How many unique features in the annotation file Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3?
19
15
17
18 

### How many lines are in the file Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3?
8936104
56976`
502932
347

### How many words are in the file Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3?  
8936104
56976
502932`
347

### How big is the file Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3 in bytes?
8936104`
56976
502932
347
