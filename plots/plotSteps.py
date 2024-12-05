import matplotlib.pyplot as plt

data = []
with open('data.txt', 'r') as f:
    for line in f:
        data.append([float(x) for x in line.split()])

data_points = []
with open('data_points.txt', 'r') as f:
    for line in f:
        data_points.append([float(x) for x in line.split()])

fig, axs = plt.subplots(4, figsize=(8, 6))
ylabels = ['x1', 'x2', 'v1', 'v2']
for i, ax in enumerate(axs):
    ax.set_ylabel(ylabels[i])

x1 = [x[0] for x in data]
x2 = [x[1] for x in data]
v1 = [x[2] for x in data]
v2 = [x[3] for x in data]

x1s = [x[0] for x in data_points]
x2s = [x[1] for x in data_points]
v1s = [x[2] for x in data_points]
v2s = [x[3] for x in data_points]
data_point_steps= [x[4] for x in data_points]

axs[0].plot(x1)
axs[0].plot(data_point_steps, x1s, 'o')
axs[1].plot(x2)
axs[1].plot(data_point_steps, x2s, 'o')
axs[2].plot(v1)
axs[2].plot(data_point_steps, v1s, 'o')
axs[3].plot(v2)
axs[3].plot(data_point_steps, v2s, 'o')

plt.show()
