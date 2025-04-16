import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

script_dir = Path(__file__).parent

csv_path = script_dir / "csv_data" / "integrate_data.csv"

data = pd.read_csv(csv_path)

x1 = data['x1'].values
v1 = data['v1'].values
x2 = data['x2'].values
v2 = data['v2'].values

fig, ax = plt.subplots()
ax.set_xlim(min(min(x1), min(x2)) - 0.1, max(max(x1), max(x2)) + 0.1)
ax.set_ylim(-0.1, 0.1)

block1, = ax.plot([], [], 'bo', markersize=15)
block2, = ax.plot([], [], 'ro', markersize=15)


def init():
    block1.set_data([], [])
    block2.set_data([], [])
    return block1, block2,


def update(frame):
    block1.set_data([x1[frame]], [0])
    block2.set_data([x2[frame]], [0])
    return block1, block2,


ani = animation.FuncAnimation(fig, update, frames=len(x1), init_func=init,
                              blit=True, interval=5)

# ani.save('animation.mp4', writer='ffmpeg', fps=60)

plt.show()
