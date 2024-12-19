import matplotlib.pyplot as plt

stateFile = 'State.txt'
diffFile =  'diffStateByParams.txt'

with open(stateFile, 'r') as f:
    raw_data = [float(line.strip()) for line in f]

with open(diffFile, 'r') as f:
    raw_data_diff = [ [float(x) for x in line.strip().split()] for line in f]

state = [[],[],[],[ ]]
for i in range(0, len(raw_data), 4):
    state[0].append(raw_data[i])
    state[1].append(raw_data[i+1])
    state[2].append(raw_data[i+2])
    state[3].append(raw_data[i+3])

stateLabels = ['v1', 'x1', 'v2', 'x2']
for i in range(4):
    plt.plot(state[i], label=f'data[{i}]')

plt.legend()
diffPlotData = [[] for i in range(4)]
for rowId in range(0, len(raw_data_diff), 4):
    for i in range(4):
        diffPlotData[i].append(raw_data_diff[rowId+i])
fig2 = plt.figure()
for stateId in range(4):
    ax = fig2.add_subplot(4,1,stateId+1)
    labels = ['v10', 'x10', 'v20', 'x20', 'k1', 'k2', 'l1', 'l2', 'm1', 'm2']
    # paramIdsToPlot = [0,  2, 4, 5, 6, 7, 8, 9]
    paramIdsToPlot = [x for x in range(10)]
    for paramId in paramIdsToPlot:
        paramData = []
        for sampleId in range(len(diffPlotData[stateId])):
            paramData.append(diffPlotData[stateId][sampleId][paramId])
        ax.plot(paramData, label=labels[paramId])
    ax.legend()
#for i in range(4):
#    plt.figure()  # Open a new plot window
#        for diffByState in diffPlotData:
 #   plt.plot(diffByState, label=f'diffByState[{i}]')

plt.show()
