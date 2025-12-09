#!/usr/bin/env ipython
import pandas as pd
import matplotlib.pyplot as plt

# data summary change range to fit dataset size
df_rank = []
for i in range(1, 5):
    dataframes = []
    for j in range(i):
        df = pd.read_csv(f"profile_{i}/Profiler{j}.csv", header=None)
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
    "CG Iteration",
    "Residual Calculation",
    "Reduction",
    "dot Product",
    "A dot Product",
]
for name, group in summary.loc[filter].groupby("name"):
    group.plot(x="Ranks", y="duration", marker="o", ax=ax, label=name)

ax.set_xlabel("#Ranks")
ax.set_ylabel("Total Duration[s]")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(bottom=0)
ax.set_title("Computation Average over Ranks[loglog]")
fig.savefig("evaluation.svg")
fig.show()


# Communication Plot
fig, ax = plt.subplots()
filter = [
    "MPI Communication Init",
    "MPI Communication Wait",
]
for name, group in summary.loc[filter].groupby("name"):
    group.plot(x="Ranks", y="duration", marker="o", ax=ax, label=name)

ax.set_xlabel("#Ranks")
ax.set_ylabel("Total Duration[s]")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(bottom=0)
ax.set_title("Communication Average over Ranks[loglog]")
fig.savefig("communication.svg")
fig.show()

# Efficiency Plot
fig, ax = plt.subplots()
runtime = summary.loc["Time Step"][summary.loc["Time Step"]["Ranks"] == 1]["duration"]
efficiency_df = summary.loc["Time Step"]
efficiency_df["efficiency"] = (
    efficiency_df["duration"] * efficiency_df["Ranks"] / runtime
)
efficiency_df.loc["Time Step"]["Ranks"] = summary.loc["Time Step"]["Ranks"]

efficiency_df.plot(
    x="Ranks", y="efficiency", marker="o", ax=ax, label="Parallel Efficiency"
)
ax.set_xlabel("#Ranks")
ax.set_ylabel("Parallel Efficiency")
# ax.set_yscale("log")
# ax.set_xscale("log")
ax.set_ylim((0.5, 1.5))
ax.set_title(r"Efficiency for a $1024^2$ Grid")
fig.savefig("efficiency.svg")
fig.show()


# group = df.groupby(by="name").agg({"calltime[ms]": "mean", "#calls": "sum"})

fig, ax = plt.subplots()
filter = [
    "Compute dt",
    "Reduction",
    "dot Product",
    "A dot Product",
    "Residual Calculation",
]
for name, group in summary.loc[filter].groupby("name"):
    group.plot(x="Ranks", y="average", marker="o", ax=ax, label=name)

ax.set_xlabel("#Ranks")
ax.set_ylabel("Average Duration[s]")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(bottom=0)
ax.set_title("Communication Average over Ranks[loglog]")
fig.savefig("reduction.svg")
fig.show()
