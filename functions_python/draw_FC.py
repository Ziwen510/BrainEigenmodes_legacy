import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
import scipy.io as sio
import os

# Robust import whether used as a module or a script
try:
    from .calc_triu_ind import calc_triu_ind
except Exception:
    from calc_triu_ind import calc_triu_ind

def load_mat_robust(filepath):
    """Robustly load MATLAB files, handling both v5 and v7.3 (HDF5) formats."""
    try:
        data = sio.loadmat(filepath)
        return data
    except NotImplementedError:
        import h5py
        with h5py.File(filepath, 'r') as f:
            data = {key: f[key][()] for key in f.keys()}
        return data

def draw_FC(data_parc, fname):
    # Load colormap
    colormap_path = os.path.join(os.path.dirname(__file__), '..', 'Rapaeh_color_table_FCD.mat')
    try:
        c_map_data = load_mat_robust(colormap_path)
        cmap_arr = c_map_data.get('c3', None)
        c_map = ListedColormap(cmap_arr) if cmap_arr is not None else plt.cm.RdBu_r
    except:
        c_map = plt.cm.RdBu_r

    data_parc = np.asarray(data_parc)
    
    # Check if input is already an FC matrix or needs correlation computation
    if data_parc.ndim == 2 and data_parc.shape[0] == data_parc.shape[1]:
        FC_emp = data_parc.astype(float)
        np.fill_diagonal(FC_emp, 1.0)
    else:
        FC_emp = np.corrcoef(data_parc)
        np.fill_diagonal(FC_emp, 1.0)

    # Clean up FC matrix
    FC_emp = np.nan_to_num(FC_emp, nan=0.0, posinf=1.0, neginf=-1.0)
    FC_emp = 0.5 * (FC_emp + FC_emp.T)

    # Compute color limits from off-diagonal values
    n = FC_emp.shape[0]
    off_diag_mask = ~np.eye(n, dtype=bool)
    off_diag_vals = FC_emp[off_diag_mask]
    max_abs = np.nanpercentile(np.abs(off_diag_vals), 99.0) if off_diag_vals.size > 0 else 1.0
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    im = ax.imshow(FC_emp, cmap=c_map, norm=norm, origin='lower', aspect='equal', interpolation='nearest')
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.title(fname)
    plt.xlabel('Parcel')
    plt.ylabel('Parcel')
    plt.axis('square')
    plt.tight_layout()

    # Save plot
    os.makedirs('outputs', exist_ok=True)
    savePath = f"./outputs/{fname}_FC.png"
    plt.savefig(savePath, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', transparent=False)
    plt.close()

    # Save upper-triangular vector
    num_parcels = FC_emp.shape[0]
    triu_ind = calc_triu_ind(np.zeros((num_parcels, num_parcels)))
    FCvec = FC_emp[triu_ind]
    sio.savemat(f"./outputs/{fname}_FC.mat", {'FCvec': FCvec})
    return FCvec