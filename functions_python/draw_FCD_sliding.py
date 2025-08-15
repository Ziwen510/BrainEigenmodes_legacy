import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from matplotlib.colors import ListedColormap
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
    except Exception as e:
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

def draw_FCD_sliding(data_parc, window_size, fname):
    try:
        # Use absolute path to the colormap file
        colormap_path = os.path.join(os.path.dirname(__file__), '..', 'Rapaeh_color_table_FCD.mat')
        c_map_data = load_mat_robust(colormap_path)
        cmap_arr = c_map_data.get('c3', None)
        c_map = ListedColormap(cmap_arr) if cmap_arr is not None else plt.cm.viridis
    except Exception as e:
        print(f"Warning: Could not load colormap file: {e}")
        c_map = plt.cm.viridis

    N, T = data_parc.shape
    tlen = T - window_size + 1
    print("FCD tlen: ")
    print(tlen)
    if tlen < 2:
        raise ValueError('window_size must be < number of timepoints.')

    nROI = data_parc.shape[0]
    fc_edge = int(nROI * (nROI - 1) / 2)
    fc_vecs = np.zeros((tlen, fc_edge))
    print("nROI: ")
    print(nROI)

    # Create mask for upper triangle (excluding diagonal)
    mask_ut = np.triu(np.ones((nROI, nROI)), k=1).astype(bool)

    for w in range(tlen):
        # Extract time window: data_parc is [nROI x T], so seg is [nROI x window_size]
        seg = data_parc[:, w:(w + window_size)]
        
        
        # We want correlations across time (window_size) for each ROI pair
        # seg is [nROI x window_size], so we want correlations across axis=1 (time)
        # np.corrcoef with rowvar=True (default) treats each row as a variable
        # Since seg is [nROI x window_size], each row is an ROI time series
        # This gives us [nROI x nROI] correlation matrix
        C = np.corrcoef(seg)
        
        # Extract upper triangle values using the mask
        fc_vecs[w, :] = C[mask_ut]

    FCD = np.corrcoef(fc_vecs)
    print("shape of Sliding FCD matrix")
    print(FCD.shape)

    plt.figure(figsize=(8, 6))
    plt.imshow(FCD, cmap=c_map, vmin=0.1, vmax=1, origin='lower')
    plt.colorbar()
    plt.title(fname)
    plt.xlabel('Window #')
    plt.ylabel('Window #')
    plt.axis('square')

    os.makedirs('./outputs', exist_ok=True)
    plt.savefig(f'./outputs/{fname}_FCD(Sliding).png', dpi=300, bbox_inches='tight')
    plt.close()
    return FCD
