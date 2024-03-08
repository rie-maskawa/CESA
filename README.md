# CESA
CESA (Coexistence-Exclusion-Synchronization-Anti-synchronization) is a Python module designed to detect significant relationships among time series sets described by integer vectors. Detailed information regarding the algorithm can be found in the associated publication.

# Requirement
* python 3.11
* poetry

# Installation
To install the package in your project created with Poetry, execute the following command:
```
poetry add git+https://github.com/rie-maskawa/CESA.git#main
```

# Usage

You can try this package using sample data provided in the example directory.

## Single Sample P-value Calculation
- First, import the modules:
```
from cesa import COEX,SYAS,combining_pval
import pandas as pd
```
<img width="299" alt="image" src="https://github.com/rie-maskawa/CESA/assets/84298724/eb6320ec-5271-47ee-91aa-49cdbd633bc1">

- Load your data file. The data should be structured as a time series matrix with species in rows and time in columns.
 ```
 data=pd.DataFrame('example/data_sample1.csv',header=0,index_col=0) 
 ```

- Use the COEX function to compute the p-values for coexistence and exclusion, as well as the measure of bias in appearance time (index V).
```
CO,EX,V=COEX(data)
```

The function returns matrices for coexitence (CO), exclusion (EX), and bias in appearance time (V). These matrices are in the form of upper triangular matrices. To access the p-value for co-occurrence between OTU2 and OTU7, for example, you would reference CO[1,6] (note that the indices may be shifted).

 - Use the SYAS function to compute the p-values for synchronization and antisynchronization, as well as the normalized inner product of the ternarized time series.
```
SY,AS,NI=SYAS(data,total_read=3000,correction=True)
```
If you analyze the relative abundance of OTUs by correcting for compositional variation, input the total read when analyzing the relative abundance of OTUs. Please input the relative abundance represented by the number of reads, not normalized abundance in the range 0-1.\
The function returns p-values for synchronization (SY), antisynchronization (AS), and the normalized inner product (NI) of the ternarized time series.

**Note**: If there is no regulation for the total sum of values, set total_read=0 and correction=False in the SYAS function.

## Combining P-values for Multiple Samples
- Create lists to store p-values and input them into the combining_pval function:
```
CO_list=[]
EX_list=[]
SY_list=[]
AS_list=[]
for s in range(S):
  data=pd.DataFrame(f'example/data_sample{s+1}.csv',header=0,index_col=0)
  CO,EX,V=COEX(data)
  SY,AS,NI=SYAS(data,total_read=3000,correction=True)
  CO_list.append(CO)
  EX_list.append(EX)
  SY_list.append(SY)
  AS_list.append(AS)
CO_all=combining_pval(CO_list)
EX_all=combining_pval(EX_list)
SY_all=combining_pval(SY_list)
AS_all=combining_pval(AS_list)
```
**Note**: Analyze common species across all samples. If there are samples with different numbers of species, it will result in an error.

# Author
* Rie Maskawa, Hideki Takayasu, Tanzila Islam, Lena Takayasu, Rina Kurokawa, Hiroaki Masuoka, Wataru Suda, Misako Takayasu
* Tokyo Institute of Technology
* masukawa.r.aa@m.titech.ac.jp

# License
CESA is under MIT license.
