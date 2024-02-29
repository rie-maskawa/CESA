# CESA
CESA (Coexistence-Exclusion-Synchronization-Anti-synchronization) is a Python module designed to detect significant relationships among time series sets described by integer vectors. Detailed information regarding the algorithm can be found in the associated publication.

# Requirement
poetry

# Installation

```
poetry add git+https://github.com/rie-maskawa/CESA.git#main
```

# Usage
```
from cesa import COEX,SYAS,combining_pval
import numpy as np
import pandas as pd

#When calculating p-values for a single sample
data=pd.DataFrame('example/data_sample1.csv',header=0,index_col=0)  #Time series matrix with species in rows and time in columns
CO,EX,ADBC=COEX(data)
SY,AS,NI=SYAS(data,total_read=3000,correction=True)  #if there is no regulation for the total sum of values, set total_read=0, correction=False.


#When calculating combined p-values for multiple samples
CO_list=[]
EX_list=[]
SY_list=[]
AS_list=[]
for s in range(S):
  data=pd.DataFrame(f'example/data_sample{s+1}.csv',header=0,index_col=0)
  CO,EX,ADBC=COEX(data)
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

# Author
* Rie Maskawa, Hideki Takayasu, Tanzila Islam, Lena Takayasu, Rina Kurokawa, Hiroaki Masuoka, Wataru Suda, Misako Takayasu
* Tokyo Institute of Technology
* masukawa.r.aa@m.titech.ac.jp

# License
"CESA" is under MIT license.
