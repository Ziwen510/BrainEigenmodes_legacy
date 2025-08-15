import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.signal import butter, filtfilt, hilbert
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

def draw_FCD_phase(ts, outName):
    try:
        # Use absolute path to the colormap file
        colormap_path = os.path.join(os.path.dirname(__file__), '..', 'Rapaeh_color_table_FCD.mat')
        c_map_data = load_mat_robust(colormap_path)
        cmap_arr = c_map_data.get('c3', None)
        c_map = plt.cm.turbo
        if cmap_arr is not None:
            try:
                cmap_arr = np.array(cmap_arr)
                cmap_arr = np.squeeze(cmap_arr)
                # Ensure shape is (N, 3)
                if cmap_arr.ndim == 2:
                    if cmap_arr.shape[1] != 3 and cmap_arr.shape[0] == 3:
                        cmap_arr = cmap_arr.T
                    # Scale if provided in 0..255
                    if np.nanmax(cmap_arr) > 1.0:
                        cmap_arr = cmap_arr / 255.0
                    if cmap_arr.shape[1] == 3:
                        c_map = ListedColormap(cmap_arr)
            except Exception:
                pass
    except Exception as e:
        print(f"Warning: Could not load colormap file: {e}")
        c_map = plt.cm.turbo

    # Match MATLAB's 1-based indexing: ts(:,22:end) -> drop first 21 columns
    ts = np.asarray(ts, dtype=float)
    ts = ts[:, 21:]
    fBand = [0.04, 0.07]
    TR = 0.72

    N, T = ts.shape

    # Match MATLAB butter design: absolute frequencies with sampling frequency fs=1/TR
    b, a = butter(2, fBand, btype='band', fs=1/TR)
    # Filter along time axis (axis=1) to mirror MATLAB's filtfilt on ts'
    tsF = filtfilt(b, a, ts, axis=1)

    # Hilbert transform along time axis (axis=1), mirroring MATLAB's hilbert(tsF.')'
    Theta = np.angle(hilbert(tsF, axis=1))

    P = int(N * (N - 1) / 2)
    Delta = np.zeros((P, T))
    idx = 0
    for i in range(N - 1):
        for j in range(i + 1, N):
            Delta[idx, :] = np.cos(Theta[i, :] - Theta[j, :])
            idx += 1

    normDelta = np.sqrt(np.sum(Delta**2, axis=0))
    FCD = (Delta.T @ Delta) / (normDelta.reshape(-1, 1) @ normDelta.reshape(1, -1))
    FCD = np.clip(FCD, -1, 1)

    plt.figure(figsize=(8, 6))
    plt.imshow(FCD, cmap=c_map, vmin=-1, vmax=1, origin='lower', interpolation='nearest')
    plt.axis('square')
    plt.colorbar()
    outName_plot = f"{outName}_FCD(Phase)"
    plt.title(outName_plot)
    plt.xlabel('Time index')
    plt.ylabel('Time index')

    os.makedirs('outputs', exist_ok=True)
    plt.savefig(f"./outputs/{outName_plot}.png", dpi=300, bbox_inches='tight')
    print(f"Saved FCD figure as '{outName_plot}'")
    plt.close()
    return FCD
