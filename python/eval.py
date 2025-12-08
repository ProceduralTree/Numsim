#!/usr/bin/env ipython
import pandas as pd
import matplotlib.pyplot as plt

df_rank = []
for i in range(1, 9):
    dataframes = []
    for j in range(i):
        df = pd.read_csv(f"profile_{i}/Profiler{j}.csv")
        df.columns = ["name", "count", "duration", "average", "min", "max"]
        for col in ["duration", "average", "min", "max"]:
            df[col] = df[col].str.replace("ns", "", regex=False).astype(float) * 1e-9
        dataframes.append(df)
    df_avg = pd.concat(dataframes).groupby(by="name").mean()
    df_avg["Ranks"] = i
    df_rank.append(df_avg)
summary = pd.concat(df_rank)
fig, ax = plt.subplots()
# Filter to include in the Plot
filter = [
    "MPI Communication Init",
    "MPI Communication Wait",
    "Time Step",
    "CG Iteration",
    "Residual Calculation",
    "dot Product",
    "Compute dt",
]
for name, group in summary.loc[filter].groupby("name"):
    group.plot(x="Ranks", y="duration", marker="o", ax=ax, label=name)

ax.set_xlabel("#Ranks")
ax.set_ylabel("Duration[s]")
ax.set_title("Duration vs Ranks per Name")
fig.savefig("evaluation.svg")
fig.show()

# group = df.groupby(by="name").agg({"calltime[ms]": "mean", "#calls": "sum"})
