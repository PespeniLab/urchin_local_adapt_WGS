# Generating allele frequencies and filtering steps

First, I run angsd on each population separately. E.g.

```bash
./angsd -b /users/c/p/cpetak/WGS/BOD_rmdups_jo.txt 
-ref /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-anc /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-out /users/c/p/cpetak/WGS/angsd_new/BOD_angsd_allsites 
-nThreads 16 
-remove_bads 1 
-C 50 
-baq 1 
-minMapQ 30 
-minQ 20 
-minInd 17 # 85% of 20 individuals
-setMinDepthInd 3 
-skipTriallelic 1 
-GL 1 
-doCounts 1 
-doMajorMinor 1 
-doMaf 1 
-doSaf 1 
-doHWE 1
```

For FOG and CAP, I  run the above code without outliers

Since population maf files were very large, I split each maf file into 1,000,000 line chuncks in separate folders inside angsd_new:

```bash
#!/bin/sh

header="chromo\tposition\tmajor\tminor\tref\tanc\tknownEM\tnInd"

while read line ; do
        pop1=$(cut -d '.' -f1 <<< $line)
        echo $pop1
        mkdir ${pop1}_split
        cp /users/c/p/cpetak/WGS/angsd_noout/${pop1}_angsd_allsites.mafs ${pop1}_split 			#angsd_new for every pop except for CAP and FOG
        cd ${pop1}_split
        head -1 ${pop1}_angsd_allsites.mafs > beg
        tail -n +2 ${pop1}_angsd_allsites.mafs > ${pop1}_tail.mafs
        split -l 1000000 ${pop1}_tail.mafs segment
        sed -i -f - segment* < <(sed 's/^/1i/' beg)
        cd ..
        sleep 0.5
done < $1
```

Once I have all the smaller files named segment in each pop_split folder I do the following:

```bash
#!/bin/sh

files=$(ls|grep segment)

for f in $files ;
do
    echo $f
    FILE=$(mktemp)
    cat /users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/header.sh >> $FILE
    echo "cd /users/c/p/cpetak/WGS/angsd_new/${1}_split" >> $FILE
    echo "python /users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/process_input.py ${f} ${1}" >> $FILE
    sbatch $FILE
    sleep 0.5
    rm $FILE
done
```

Where process_input.py:

```python
import pandas as pd
import os
import numpy as np
import sys

file_name = sys.argv[1]
pop = sys.argv[2]

file = '/users/c/p/cpetak/WGS/angsd_new/' + pop + '_split/' + file_name

data = pd.read_csv(file, sep="\t")
data['filename'] = pop

print("read data")

def init_proc_file(data):
    data['chromo'] = data['chromo'].apply(str)
    data['position'] = data['position'].apply(str)
    data['pos'] = data[['chromo', 'position']].apply(lambda x: ''.join(x), axis=1)
    data.loc[data.minor == "A", 'A_freq'] = data.knownEM * data.nInd * 2
    data.loc[data.minor == "T", 'T_freq'] = data.knownEM * data.nInd * 2
    data.loc[data.minor == "C", 'C_freq'] = data.knownEM * data.nInd * 2
    data.loc[data.minor == "G", 'G_freq'] = data.knownEM * data.nInd * 2
    data.loc[data.major == "A", 'A_freq'] = (data.nInd*2)-(data.knownEM * data.nInd * 2)
    data.loc[data.major == "T", 'T_freq'] = (data.nInd*2)-(data.knownEM * data.nInd * 2)
    data.loc[data.major == "C", 'C_freq'] = (data.nInd*2)-(data.knownEM * data.nInd * 2)
    data.loc[data.major == "G", 'G_freq'] = (data.nInd*2)-(data.knownEM * data.nInd * 2)
    data = data.drop(['major'], axis=1)
    data = data.drop(['minor'], axis=1)
    data = data.drop(['ref'], axis=1)
    data = data.drop(['anc'], axis=1)
    data = data.drop(['nInd'], axis=1)
    data = data.drop(['chromo'], axis=1)
    data = data.drop(['position'], axis=1)
    data = data.drop(['knownEM'], axis=1)
    pop = data['filename'].iloc[1]
    print(pop)
    data = data.drop(['filename'], axis=1)
    data = data.rename(columns={'A_freq': 'A_freq'+pop, 'T_freq': 'T_freq'+pop,'C_freq': 'C_freq'+pop,'G_freq': 'G_freq'+pop})
    data = data.round()
    return(data)


new_df = init_proc_file(data)

outfile_name = '/users/c/p/cpetak/WGS/angsd_new/' + pop + '_split/' + file_name + '_processed.csv'
new_df.to_csv(outfile_name, index=False)
```

to get a processed file for each segment. then, I merge these files into one:

```bash
cat *processed.csv > combined.csv
grep -v 'pos' combined.csv > combined_fixed.csv
head -1 segmentaa_processed.csv > header.txt
cat header.txt combined_fixed.csv >> combined_fixed2.csv

cd my_bayenv/testing_new
cp ../../TER_split/combined_fixed2.csv ./TER # for all pops
```

Then merged all populations into 1 using the following code (merge_dd.py):

```python
import pandas as pd
import os
import numpy as np
import dask.dataframe as dd
import sys

file1 = '/users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/' + sys.argv[1]
file2 = '/users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/' + sys.argv[2]

data1 = dd.read_csv(file1, sep=",")
data2 = dd.read_csv(file2, sep=",")

print("read data")

merged = dd.merge(data1, data2, on="pos", how="inner")

print("merged")

outfile_name = '/users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/' + sys.argv[3]

merged.to_csv(outfile_name, index=False)
```

NOTE: this way only sites where every pop has information about it is kept!

NOTE: since the output was a directory not a csv, I had to do the following each time:

```bash

cat * > combined
cp 1.csv/combined ./1_c.csv
rm -r 1.csv
head -1 1_c.csv > fline
grep -v 'pos' 1_c.csv > 1.csv
cat fline 1.csv > 1_new.csv
rm 1.csv
rm 1_c.csv
rm fline
mv 1_new.csv 1.csv
```

then, I processed this merged file (in separate chunks as it was too big):

```python
import pandas as pd
import numpy as np
import sys

print(sys.argv[0])
file_name = sys.argv[1]

merged=pd.read_csv(file_name)
merged2=merged.fillna(0)

def flatten(l):
    return [i for seq in l for i in seq]

def divide_chunks(l, chunk_size):
    # looping till length l
    for i in range(0, len(l), chunk_size):
        yield l[i : i + chunk_size]

def unzip(l):
    # turn list of pairs into pairs of lists
    return list(zip(*l))

npops = 7
column_names = ["pos", "LOMal1","LOMal2","SANal1","SANal2","TERal1","TERal2","BODal1","BODal2","KIBal1","KIBal2","CAPal1","CAPal2","FOGal1","FOGal2"]
new_df=pd.DataFrame(columns = column_names)

for i, row in merged2.iterrows():
    values = list(row)
    popinfo = values[1 : 1 + npops * 4]
    otherinfo = values[0:1]
    notnans = [str(val) != "0.0" for val in popinfo]
    notnansperpop = list(divide_chunks(notnans, 4))
    #print(notnansperpop)
    colstocompare = unzip(notnansperpop)
    #print(colstocompare)
    nonzeros=[True in set(vals) for vals in colstocompare]
    #print(nonzeros)
    num_nonzero=nonzeros.count(True)
    #print(num_nonzero)
    if num_nonzero == 2:
      mask = nonzeros
      #print(mask)
      pops = list(divide_chunks(popinfo, 4))
      #print(pops)
      final = flatten([np.array(p)[mask] for p in pops])
      #print(final)
      list_to_add=otherinfo+final
      #print(list_to_add)
      new_df.loc[len(new_df)] = list_to_add
      #print("success")

outfile_name=file_name + "_processed.csv"

new_df.to_csv(outfile_name, index=False)
```

output after combining split files: merged_processed.csv

Filtering for MAF (process_final_withfilter.py)

```python
import pandas as pd
import numpy as np

merged=pd.read_csv("merged_processed.csv")

MAF=0.025
MAC=MAF*274 #because we dropped 3 individuals so 6 chromosomes

columns = ['LOMal1', 'SANal1', 'TERal1', 'BODal1','KIBal1','CAPal1','FOGal1']
merged2 = merged[merged[columns].sum(axis=1) > MAC]
columns2 = ['LOMal2', 'SANal2', 'TERal2', 'BODal2','KIBal2','CAPal2','FOGal2']
merged3 = merged2[merged2[columns2].sum(axis=1) > MAC]


final_df = pd.DataFrame([el for el in merged3["pos"].values for _ in range(2)])
#final_df = final_df.rename(columns={0: "pos"})
final_df["LOM"] = [ el for pair in zip(merged3["LOMal1"].values, merged3["LOMal2"].values) for el in pair]
final_df["SAN"] = [ el for pair in zip(merged3["SANal1"].values, merged3["SANal2"].values) for el in pair]
final_df["TER"] = [ el for pair in zip(merged3["TERal1"].values, merged3["TERal2"].values) for el in pair]
final_df["BOD"] = [ el for pair in zip(merged3["BODal1"].values, merged3["BODal2"].values) for el in pair]
final_df["KIB"] = [ el for pair in zip(merged3["KIBal1"].values, merged3["KIBal2"].values) for el in pair]
final_df["CAP"] = [ el for pair in zip(merged3["CAPal1"].values, merged3["CAPal2"].values) for el in pair]
final_df["FOG"] = [ el for pair in zip(merged3["FOGal1"].values, merged3["FOGal2"].values) for el in pair]

final_df.to_csv('/users/c/p/cpetak/WGS/angsd_new/my_bayenv/testing_new/bayenv_withpos_0025filter.csv', sep='\t', index=False)
```

Final number of sites: 994,220

No outliers, biallelic (across all pops, not just per pop), minMAF 0.025, every pop has site information





