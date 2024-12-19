import matplotlib.pyplot as plt

fileName = 'resid.txt'
with open(fileName, 'r') as f:
    raw_data = [float(line.strip()) for line in f]

state = [[],[],[],[ ]]
for i in range(0, len(raw_data), 4):
    state[0].append(raw_data[i])
    state[1].append(raw_data[i+1])
    state[2].append(raw_data[i+2])
    state[3].append(raw_data[i+3])

stateLabels = ['v1', 'x1', 'v2', 'x2']
for i in range(4):
    plt.plot(state[i], label=stateLabels[i])

plt.legend()
plt.show()