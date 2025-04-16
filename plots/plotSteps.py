import matplotlib.pyplot as plt

def loadStateData(fileName):
    
    with open(fileName, 'r') as f:
        raw_data = [float(line.strip()) for line in f]

    state = [[],[],[],[ ]]
    for i in range(0, len(raw_data), 4):
        state[0].append(raw_data[i])
        state[1].append(raw_data[i+1])
        state[2].append(raw_data[i+2])
        state[3].append(raw_data[i+3])
    return state

if __name__ == "__main__":
    state = loadStateData('../build/data_points.txt')
    optState = loadStateData('../build/opt_state_data_points.txt')
    stateLabels = ['v1', 'x1', 'v2', 'x2']
    fig, axs = plt.subplots(4, figsize=(8, 6))
    for i in range(4):
        axs[i].plot(state[i], label=stateLabels[i])
        axs[i].plot(optState[i], label=stateLabels[i])
        axs[i].set_ylabel(stateLabels[i])
    
    plt.legend()
    plt.show()