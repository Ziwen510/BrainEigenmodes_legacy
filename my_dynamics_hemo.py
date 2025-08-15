import os
from pathlib import Path
import numpy as np
import scipy.io as sio

try:
    import h5py  # optional for v7.3 MAT support
except Exception:
    h5py = None

from analysis_pipeline import analysis_pipeline


def load_tc_mat(filename: str) -> np.ndarray:
    """
    Load the time-course matrix 'tc' from a MAT-file and return as [vertices x time].
    Tries scipy.io.loadmat first; on failure or v7.3 files, falls back to h5py.
    """
    filename = str(filename)
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    # Try standard MAT load
    try:
        S = sio.loadmat(filename)
        if 'tc' not in S:
            raise KeyError("Variable 'tc' not found in MAT file")
        tc = S['tc']  # typically [T x N]
        if tc.ndim != 2:
            tc = np.squeeze(tc)
            if tc.ndim != 2:
                raise ValueError("Loaded 'tc' does not have 2 dimensions")
        return tc.T  # return [N x T]
    except Exception:
        # Fall back to HDF5
        if h5py is None:
            raise
        with h5py.File(filename, 'r') as f:
            if 'tc' in f:
                tc = f['tc'][...]
            elif '/tc' in f:
                tc = f['/tc'][...]
            else:
                # guess first dataset
                key = next(iter(f.keys()))
                tc = f[key][...]
        return np.array(tc).T


def main():
    # Parameters mirroring the MATLAB script
    surface_interest = 'fsLR_32k'
    hemisphere = 'lh'

    # Load cortex mask
    cortex_mask_path = Path('BrainEigenmodes_legacy/data/template_surfaces_volumes') / f"{surface_interest}_cortex-{hemisphere}_mask.txt"
    if not cortex_mask_path.exists():
        # try relative to current script directory
        here = Path(__file__).resolve().parent
        cortex_mask_path = here / 'data/template_surfaces_volumes' / f"{surface_interest}_cortex-{hemisphere}_mask.txt"
    cortex = np.loadtxt(cortex_mask_path).astype(bool)
    cortex_ind = np.where(cortex)[0]
    print(f"Number of vertices in cortex: {cortex_ind.size}")

    # Timing parameters (ms)
    tmax = 865 * 1000
    tstep = 720
    print("specified Nsteps:", int(tmax / tstep))

    # Naming and model params
    input_mode = 'noise'  # or 'pulse'
    Tmax = tmax / 1000  # seconds
    simulator = 'cortical'
    norm = ''
    version = '_fused_acc'
    rmax = '430'
    variable = '_bold'
    boldnorm = ''
    cropping = ''

    rate = 1000
    nu_es = 0.0001 * rate
    nu_ee = 0.00006 * rate
    v0 = 0.0006 * rate
    theta = 0.0126766 * rate
    sigma = 0.0038 * rate
    r_s = 18
    v = r_s * 116

    if input_mode == 'noise':
        seed = 1
        mean_I = 0.0
        ASD = 0.01
        str_part = f"white_ASD{ASD:.4f}_seed{seed}_mean{mean_I:.1f}"
    else:  # 'pulse'
        cite = 13289
        amp = 1000
        str_part = f"pulse_{cite}_I{amp}"

    # Build base filename consistent with simulator outputs
    base_fname = (
        f"cortical{version}_{str_part}_"
        f"nuee{nu_ee:.4f}_nues{nu_es:.4f}_r{r_s}_v{v}_v0{v0:.4f}_"
        f"theta{theta:.4f}_sigma{sigma:.4f}_T{Tmax:.0f}_rmax{rmax}{norm}{boldnorm}"
    )

    # Mat path (load from NFT/mats)
    mats_dir = Path('/home/zwang/storage/my/NFT/mats')
    mat_path = mats_dir / f"{base_fname}{variable}.mat"

    print(f"Loading BOLD MAT: {mat_path}")
    simulated_activity_rest = load_tc_mat(mat_path)
    print("Size of loaded TC:", simulated_activity_rest.shape)

    # Downsample/select time window
    warmup = 1440  # ms
    target_tstep = 720  # ms

    if cropping == "_cropped":
        T_target = np.arange(warmup, 217*1000 + 1, target_tstep)
    else:
        T_target = np.arange(warmup, tmax + 1, target_tstep)

    idx = (np.round(T_target / tstep)).astype(int)  # MATLAB had +1 for 1-based indices
    tc_down = simulated_activity_rest[:, idx]

    # Embed into full vertex array
    nc, Tn = tc_down.shape
    Ntotal = cortex.size
    full_activity = np.zeros((Ntotal, Tn), dtype=tc_down.dtype)
    full_activity[cortex_ind, :] = tc_down
    print("size of final full tc:", full_activity.shape)

    # Build output name
    out_name = f"my_{simulator}{version}_{str_part}_r{r_s}_T{Tmax:.1f}{norm}{cropping}_{variable}"
    print("Output name:", out_name)

    # Call analysis pipeline (convert ms to seconds for plotting)
    T_seconds = T_target / 1000.0
    metrics = analysis_pipeline(full_activity, out_name, T_seconds, cortex)
    print(metrics)


if __name__ == "__main__":
    main()
