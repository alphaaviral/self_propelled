import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("./out.csv", header=None, index_col=0)
# drop index column
# df.drop(df.columns[0], axis=1, inplace=True)
arr=df.to_numpy()
arr=arr.reshape([5,22,6])

for i in range(0,5):
    plt.clf()
    for j in range(0,22):
        plt.plot(arr[i][j][0], arr[i][j][1], 'r', marker=".", markersize=6)
    plt.xlim(-50, 50)
    plt.ylim(-50, 50)
    plt.savefig("./test"+str(i)+".png")
    print(i)

# print(df)