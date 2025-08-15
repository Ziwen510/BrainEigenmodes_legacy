import numpy as np
import matplotlib.pyplot as plt
import os

def draw_tc_samples(simulated_activity_rest, cortex_ind, name, T, num_samples=6, is_cropped=False):
    cortex_ind = np.array(cortex_ind).flatten()
    cortex_ind = cortex_ind[(cortex_ind >= 1) & (cortex_ind <= simulated_activity_rest.shape[0])]

    Nc = len(cortex_ind)
    if num_samples > Nc:
        print(f'Warning: Requested {num_samples} samples but only {Nc} cortex vertices available. Reducing to {Nc}.')
        num_samples = Nc

    sample_positions = np.round(np.linspace(1, min(20000, Nc), num_samples)).astype(int) - 1
    idx_samples = cortex_ind[sample_positions]

    n = len(idx_samples)
    n_cols = 2
    n_rows = int(np.ceil(n / n_cols))

    if is_cropped:
        T = T[300:]
        simulated_activity_rest = simulated_activity_rest[:, 300:]
        name = f"{name}_cropped"

    fig_raw = plt.figure(figsize=(12, 4*n_rows))
    fig_raw.suptitle(name, fontsize=16)

    for k in range(n):
        vi = idx_samples[k]
        plt.subplot(n_rows, n_cols, k + 1)
        plt.plot(T, simulated_activity_rest[vi-1, :], linewidth=1)
        plt.xlim([T[0], T[-1]])
        plt.xlabel('Time (s)')
        plt.ylabel('Activity')
        plt.title(f'Vertex {vi}')
        plt.grid(True)

    plt.tight_layout()
    os.makedirs('outputs', exist_ok=True)
    plt.savefig(os.path.join('outputs', f'{name}.png'), dpi=300, bbox_inches='tight')
    plt.close()
