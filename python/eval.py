#!/usr/bin/env ipython
import pandas as pd
import matplotlib.pyplot as plt

sizes = [256, 512, 1024]
# sizes = [64, 128, 256]

# data summary change range to fit dataset size
summaries = []
for n in sizes:
    df_rank = []
    for i in range(1, 9):
        dataframes = []
        for j in range(i):
            df = pd.read_csv(f"profile_mpi{i}_n{n}_{n}/Profiler{j}.csv", header=None)
            df.columns = ["name", "count", "duration", "average", "min", "max"]
            for col in ["duration", "average", "min", "max"]:
                df[col] = (
                    df[col].str.replace("ns", "", regex=False).astype(float) * 1e-9
                )

            dataframes.append(df)

        df_avg = pd.concat(dataframes).groupby(by="name").mean()
        df_avg["Ranks"] = i
        df_rank.append(df_avg)
    summary = pd.concat(df_rank)
    summary["nCells"] = n
    summaries.append(summary)


for summary in summaries:
    for name, group in summary.groupby("Ranks"):
        runtime = summary[summary["Ranks"] == 1]
        summary.loc[group.index, "average efficiency"] = (
            runtime["average"] / summary["average"]
        )
        summary.loc[group.index, "efficiency"] = (
            runtime["duration"] / summary["duration"]
        )


for summary, n in zip(summaries, sizes):
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
        group.plot(x="Ranks", y="average efficiency", marker="o", ax=ax, label=name)

    ax.set_xlabel("#Ranks")
    ax.set_ylabel("Average Speedup")
    ax.set_ylim(bottom=0)
    ax.set_title(f"Computation Speedup for a {n}x{n} Grid")
    fig.savefig(f"images/computation_speedup_n{n}_{n}.svg")

for summary, n in zip(summaries, sizes):
    fig, ax = plt.subplots()

    # Filter to include in the Plot
    filter = [
        "MPI Communication Init",
        "MPI Communication Wait",
    ]
    for name, group in summary.loc[filter].groupby("name"):
        group.plot(x="Ranks", y="average efficiency", marker="o", ax=ax, label=name)

    ax.set_xlabel("#Ranks")
    ax.set_ylabel("Average Speedup")
    ax.set_ylim(bottom=0)
    ax.set_title(f"Communication Speedup for a {n}x{n} Grid")
    fig.savefig(f"images/communication_speedup_n{n}_{n}.svg")

for summary, n in zip(summaries, sizes):
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
    ax.set_title(f"Computation Average for a {n}x{n} Grid")
    fig.savefig(f"images/evaluation_n{n}_{n}.svg")

for summary, n in zip(summaries, sizes):
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
    ax.set_ylim(bottom=0)
    ax.set_title(f"Communication Average for a {n}x{n} Grid")
    fig.savefig(f"images/communication_n{n}_{n}.svg")


for summary, n in zip(summaries, sizes):
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
    ax.set_title(f"Communication Average for a {n}x{n} Grid")
    fig.savefig(f"images/reduction_n{n}_{n}.svg")


# Efficiency Plot
fig, ax = plt.subplots()
for summary, n in zip(summaries, sizes):
    runtime = summary.loc["Time Step"][summary.loc["Time Step"]["Ranks"] == 1][
        "duration"
    ]
    efficiency_df = summary.loc["Time Step"]
    efficiency_df["efficiency"] = runtime / (
        efficiency_df["duration"] * efficiency_df["Ranks"]
    )
    efficiency_df.loc["Time Step"]["Ranks"] = summary.loc["Time Step"]["Ranks"]

    efficiency_df.plot(
        x="Ranks", y="efficiency", marker="o", ax=ax, label=f"Gridsize {n}x{n}"
    )
ax.set_xlabel("#Ranks")
ax.set_ylabel("Parallel Efficiency")
ax.set_ylim((0.0, 2))
ax.set_title(f"Parallel Efficiency")
fig.savefig("images/efficiency.svg")
