#### Which one of the following commands is used for replacing "|" with ":" globally in test file1?
sed 's/|/:/g' file1 `
sed '/|/:/g' file1
sed 's/\|/\:/' file1
sed 's/|/:/' file1

#### There is a sample file, fruit.txt, as shown below:

### apple

### orange

### banana

### papaya

#### Which one of the following commands is used to add "Fruit:" to the beginning of each line in the fruit.txt file?
sed 's/^/Fruit:/' fruit.txt `
sed 's/$/Fruit:/' fruit.txt
sed 's/^Fruit://' fruit.txt
sed 's/$Fruit://' fruit.txt

#### There is a sample file, fruit.txt, as shown below:

### apple

### orange

### banana

### papaya

#### Which one of the following commands is used to add "sweet" to the end of each line in the fruit.txt file?
sed 's/^/sweet/' fruit.txt 
sed 's/$/sweet/' fruit.txt `
sed 's/^sweet//' fruit.txt
sed 's/$sweet//' fruit.txt



#### There is a sample file, employee.txt, shown as below:

Ajay manager account 45000

Sunil clerk account 25000

Varun manager sales 50000

Amit manager account 47000

Tarun assistant sales 15000

Deepak clerk sales 23000

Satvik director purchasing 80000 

Sunil assistant  sales 13000

#### Which one of the following commands is used to print the employee name and salary fields?
awk '{print $1&&$4}' employee.txt
awk '{print $0,$1,$4}' employee.txt
awk '{print $1&$4}' employee.txt
awk '{print $1,$4}' employee.txt


#### Which one of the following commands can not print the salary less than 50000 for the employee whose title is manager?
awk '$2 == "manager" && $4 < 50000' employee.txt
awk '/manager/ && $4 < 50000' employee.txt
awk '$2 = "manager" && $4 < 50000' employee.txt `
awk '$2 ~ /manager/ && $4 < 50000' employee.txt


#### Which one of the following commands can print the employee name whose salary is more than 60000?
awk '$4 > 60000 && print $1' employee.txt
awk '$4 > 60000' '{print $1}' employee.txt
awk '$4 > 60000 {print $1}' employee.txt `	
awk '$4 > 60000 && $1' employee.txt


#### For the file gff3file (Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3), how many mRNAs have a starting position after 20,000,000 bases and are on the + strand? 
787`
1691
52036
55091

#### For the file gff3file (Monodelphis_domestica.ASM229v1.101.primary_assembly.7.gff3), how many exons have a length more than 2000 nt?
198`
4890
1523
29875

#### Which one of the following commands can print out the variable my_var?
echo my_var
echo $my_var
cat my_var
less my_var


#### Can you think of where we could put the > so that it wouldn't overwrite the file with each iteration of the loop? 
for item in {car,truck,ukulele}; echo $item > words.txt; done
for item in {car,truck,ukulele}; do echo >> $item ; done
for item in {car,truck,ukulele}; echo $item >> words.txt; done 
for item in {car,truck,ukulele}; do echo $item; done > words.txt` 


