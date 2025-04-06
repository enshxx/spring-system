import pandas as pd
import matplotlib.pyplot as plt


def show_plt(csv_path):
    data = pd.read_csv(csv_path)

    x1 = data['x1'].values
    v1 = data['v1'].values
    x2 = data['x2'].values
    v2 = data['v2'].values

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(min(min(x1), min(x2)) - 0.5, max(max(x1), max(x2)) + 0.5)
    ax.set_ylim(min(min(v1), min(v2)) - 0.5, max(max(v1), max(v2)) + 0.5)

    block1, = ax.plot(x1, v1, 'bo', markersize=5, label='Грузик 1')
    block2, = ax.plot(x2, v2, 'ro', markersize=5, label='Грузик 2')

    ax.set_xlabel('Положение, x', fontsize=12)
    ax.set_ylabel('Скорость, v', fontsize=12)
    ax.set_title('График зависимости скорости от положения', fontsize=14, pad=20)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(fontsize=12)

    plt.show()
