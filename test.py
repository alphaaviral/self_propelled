import pandas as pd
import numpy as np
import csv
test = np.arange(12000000).reshape([10000,200,6])
# test = [[[[1,2,3],[4,5,6]],
        # [[7,8,9],[10,11,12]],
        # [[13,14,15],[16,17,18]]],

        # [[[19,20,21],[22,23,24]],
        # [[25,26,27],[28,29,30]],
        # [[31,32,33],[34,35,36]]]]
# test=np.array(test)
newtest=test.reshape([10000,1200])
# newtest=newtest.reshape([10000,200,6])
df=pd.DataFrame(newtest)
df.to_csv("./test.csv",header=False)
names = ['t','chain', 'particle', 'attribute']
index = pd.MultiIndex.from_product([range(s)for s in test.shape], names=names)
print(index)
df = pd.DataFrame({'test': test.flatten()}, index=index)['test']
df = df.unstack(level='attribute').sort_index()
df.columns = ['A', 'B', 'C','D','E','F']
# df.index.names = ['DATE', 'i']
print(df)
df.to_csv("./test.csv")