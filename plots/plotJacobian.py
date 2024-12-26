import matplotlib.pyplot as plt

diffFile =  'build/jacobian.txt'

with open(diffFile, 'r') as f:
    raw_data_diff = [ [float(x) for x in line.strip().split()] for line in f]

diffPlotData = [[] for i in range(4)]
for rowId in range(0, len(raw_data_diff), 4):
    for i in range(4):
        diffPlotData[i].append(raw_data_diff[rowId+i])

fig = plt.figure()
for stateId in range(4):
    ax = fig.add_subplot(4,1,stateId+1)
    labels = ['v10', 'x10', 'v20', 'x20', 'k1', 'k2', 'l1', 'l2', 'm1', 'm2']
    for paramId in range(10):
        paramData = []
        for sampleId in range(len(diffPlotData[stateId])):
            paramData.append(diffPlotData[stateId][sampleId][paramId])
        ax.plot(paramData, label=labels[paramId])
    ax.legend()

plt.show()
