import matplotlib.pyplot as plt

for paramId in range(10):

    plus_data = []
    with open('plus_data' + str(paramId) + '.txt', 'r') as f:
        for line in f:
            plus_data.append([float(x) for x in line.split()])

    minus_data = []
    with open('minus_data' + str(paramId) + '.txt', 'r') as f:
        for line in f:
            minus_data.append([float(x) for x in line.split()])

    fig, axs = plt.subplots(4, figsize=(8, 6))
    ylabels = ['x1', 'x2', 'v1', 'v2']
    for i, ax in enumerate(axs):
        ax.set_ylabel(ylabels[i])

    x1 = [x[0] for x in plus_data]
    x2 = [x[1] for x in plus_data]
    v1 = [x[2] for x in plus_data]
    v2 = [x[3] for x in plus_data]

    x1s = [x[0] for x in minus_data]
    x2s = [x[1] for x in minus_data]
    v1s = [x[2] for x in minus_data]
    v2s = [x[3] for x in minus_data]


    axs[0].plot(x1, 'r')
    axs[0].plot(x1s, 'b')
    axs[1].plot(x2, 'r')
    axs[1].plot(x2s, 'b')
    axs[2].plot(v1, 'r')
    axs[2].plot(v1s, 'b')
    axs[3].plot(v2, 'r')
    axs[3].plot(v2s, 'b')

plt.show()
