import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os

def load_mat_robust(filepath):
    """
    Robustly load MATLAB files, handling both v5 and v7.3 (HDF5) formats.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    try:
        # Try standard scipy.io.loadmat first (for v5 files)
        data = sio.loadmat(filepath)
        print(f"Loaded {filepath} using scipy.io.loadmat")
        return data
    except NotImplementedError as e:
        if "v7.3" in str(e):
            print(f"Standard load failed: {e}\nFalling back to h5py...")
            try:
                import h5py
                with h5py.File(filepath, 'r') as f:
                    # Convert h5py dataset to numpy array
                    data = {}
                    for key in f.keys():
                        data[key] = f[key][()]
                print(f"Loaded {filepath} using h5py")
                return data
            except Exception as h5e:
                raise IOError(f"HDF5 read failed for {filepath}:\n{h5e}")
        else:
            raise e
    except Exception as e:
        raise IOError(f"Failed to load {filepath}: {e}")

def draw_FC_corr_scatter(sim_FC, fname, emp_version):
    if emp_version == "ICA":
        s2 = load_mat_robust('/home/zwang/storage/my/BrainEigenmodes_legacy/data/empirical/group_FC_ICA.mat')
        y = s2['group_fc_mat']
    elif emp_version == "ICA+GSR":
        s2 = load_mat_robust('/home/zwang/storage/my/BrainEigenmodes_legacy/data/empirical/group_FC_ICA+GSR.mat')
        y = s2['group_fc_mat']
    elif emp_version == "ICA+GSR+WM_CSF_MT_CEN":
        s2 = load_mat_robust('/home/zwang/storage/my/BrainEigenmodes_legacy/data/empirical/group_FC_ICA+GSR+WM_CSF_MT_CEN.mat')
        y = s2['group_fc_mat']
    elif emp_version == "Minimal":
        s2 = load_mat_robust('/home/zwang/storage/my/BrainEigenmodes_legacy/data/results/model_results_Glasser360_lh.mat')
        y = s2['FC_emp']

    # Robust import whether used as a module or a script
    try:
        from .calc_triu_ind import calc_triu_ind
    except Exception:
        from calc_triu_ind import calc_triu_ind

    if sim_FC.ndim == 2 and sim_FC.shape[0] > 1:
        sim_FC = sim_FC[calc_triu_ind(sim_FC)]
    if y.ndim == 2 and y.shape[0] > 1:
        y = y[calc_triu_ind(y)]

    x = np.arctanh(sim_FC)
    y = np.arctanh(y)

    R = np.corrcoef(x, y)
    corrValue = R[0, 1]

    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.6)

    p = np.polyfit(x, y, 1)
    x_sorted = np.sort(x)
    y_fit = np.polyval(p, x_sorted)
    plt.plot(x_sorted, y_fit, 'r-', linewidth=2)

    plt.xlabel("Simulated FC")
    plt.ylabel("Empirical FC")
    plt.grid(True)
    plt.title(f'{fname} vs. {emp_version}')
    plt.text(0.05, 0.95, f'r = {corrValue:.2f}', transform=plt.gca().transAxes, verticalalignment='top', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', edgecolor='black'))

    os.makedirs('outputs', exist_ok=True)
    plt.savefig(f'./outputs/{fname}_FC_corr({emp_version}).png', dpi=300, bbox_inches='tight')
    plt.close()
    return corrValue
