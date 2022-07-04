import json
import matplotlib.pyplot as plt
import numpy
import pdb


results = {}
results2 = {}
fullResults = json.loads(
    open("Full_Benchmark_Results.json", "r").read())


results["ortools"] = json.loads(
    open("ortoolsResults.json", "r").read())
results["localsolver"] = json.loads(
    open("localsolverResults.json", "r").read())
# results["ortools"] = json.loads(
#     open("RUN_2022-03-14_local_localsolver.json", "r").read())
# results["localsolver"] = json.loads(
#     open("RUN_2022-03-14_local_ortools.json", "r").read())
REPEAT = 0
colors = {"ortools": "lime", "localsolver": "red"}
stats = ["unassigned", "cost", "total_time", "total_distance"]
# , "cost",
#  "total_time", "total_distance"]


fig = plt.figure()

axs = fig.subplots(3, 2, sharex='col')
axs = axs.flat

for statistic in stats:
    data_y = []
    data_x = None
    ax = next(axs)
    gap = []
    for solver, result in results.items():
        data_x = result["stats"].keys()
        data_y.append(
            # * (1.5 if solver=="localsolver" else 1)
            [list(stat.values())[REPEAT][statistic]
             for stat in result["stats"].values()]
        )
        ax.plot(
            data_x,
            data_y[-1],
            linestyle='none',
            marker='o',
            c=colors[solver]
        )
    if len(data_y) == 2:
        ax.fill_between(data_x,
                        [q if q != None else 0 for q in data_y[0]],
                        [q if q != None else 0 for q in data_y[1]], color='grey')
    ax.set(ylabel=statistic)


fig.tight_layout()
plt.show()


# for opt in obj_python["options"]:
#     print(opt)

# for stat in obj_python["stats"]:
#     listInstances.append(stat)
#     print(obj_python["stats"][stat])
