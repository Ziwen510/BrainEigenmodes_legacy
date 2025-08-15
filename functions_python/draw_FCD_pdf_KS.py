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
        # squeeze_me simplifies shapes like (1, N) -> (N,)
        data = sio.loadmat(filepath, squeeze_me=True, struct_as_record=False)
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

def draw_FCD_pdf_KS(mat, version, name):
    binEdges = np.linspace(0, 1, 101)
    binCenters = (binEdges[1:] + binEdges[:-1]) / 2

    utValues = mat[np.triu_indices_from(mat, k=1)]
    pdf1, _ = np.histogram(utValues, binEdges, density=True)

    # Resolve empirical PDF path robustly (prefer relative to this file)
    base_dir = os.path.dirname(__file__)
    filename = f"group_FCD_pdf_{version}.mat"
    candidate_paths = [
        os.path.normpath(os.path.join(base_dir, '..', 'data', 'empirical', filename)),
        os.path.normpath(os.path.join(os.getcwd(), 'data', 'empirical', filename)),
        os.path.normpath(os.path.join(os.getcwd(), 'BrainEigenmodes_legacy', 'data', 'empirical', filename)),
    ]
    otherPdfMatFile = None
    for p in candidate_paths:
        if os.path.isfile(p):
            otherPdfMatFile = p
            break
    if otherPdfMatFile is None:
        raise FileNotFoundError(f"Empirical PDF MAT not found. Tried: {candidate_paths}")

    S = load_mat_robust(otherPdfMatFile)

    # Try to robustly extract the empirical PDF array from the loaded dict
    target_len = len(binCenters)
    pdf2 = None

    # 1) Direct 'pdf' key
    if isinstance(S, dict) and 'pdf' in S and isinstance(S['pdf'], np.ndarray):
        pdf2 = np.asarray(S['pdf']).ravel()
    else:
        # 2) Prefer keys that look like PDF/density (exclude x-centers/edges)
        preferred_candidates = []
        fallback_candidates = []
        if isinstance(S, dict):
            for k, v in S.items():
                if not isinstance(k, str) or k.startswith('__'):
                    continue
                if not (isinstance(v, np.ndarray) and v.dtype.kind in 'iuf'):
                    continue
                arr = np.asarray(v).ravel()
                if arr.size != target_len:
                    continue
                k_lower = k.lower()
                if any(tag in k_lower for tag in ['center', 'centers', 'edge', 'edges', 'bin', 'bins', 'x_center', 'xcent']):
                    continue
                if any(tag in k_lower for tag in ['pdf', 'dens', 'prob', 'hist']):
                    preferred_candidates.append((k, arr))
                else:
                    fallback_candidates.append((k, arr))

        chosen = None
        if preferred_candidates:
            chosen = preferred_candidates[0]
        elif fallback_candidates:
            chosen = fallback_candidates[0]

        if chosen is not None:
            k_best, arr_best = chosen
            pdf2 = arr_best
            print(f"Using key '{k_best}' from {otherPdfMatFile} as empirical PDF (size={arr_best.size}).")

    if pdf2 is None:
        # 3) As a last resort, try to coerce any array-like to numeric
        for k, v in (S.items() if isinstance(S, dict) else []):
            if isinstance(k, str) and k.startswith('__'):
                continue
            try:
                arr = np.array(v).astype(float).ravel()
                if arr.size >= 10:
                    pdf2 = arr
                    print(f"Coerced key '{k}' to numeric empirical PDF (size={arr.size}).")
                    break
            except Exception:
                continue

    if pdf2 is None:
        raise ValueError(
            f"Could not locate a numeric empirical PDF in {otherPdfMatFile}. Keys: "
            f"{list(S.keys()) if isinstance(S, dict) else type(S)}"
        )

    if len(pdf2) != len(pdf1):
        raise ValueError(f'Loaded PDF length ({len(pdf2)}) does not match expected number of bins ({len(pdf1)}).')

    binWidth = binEdges[1] - binEdges[0]
    # Heuristic: if sum close to 1, treat as probability per bin and convert to density
    s = float(np.sum(pdf2))
    if 0.9 <= s <= 1.1:
        pdf2 = pdf2 / binWidth
    cdf1 = np.cumsum(pdf1) * binWidth
    cdf2 = np.cumsum(pdf2) * binWidth

    ksStat = np.max(np.abs(cdf1 - cdf2))

    plt.figure(figsize=(8, 6))
    plt.plot(binCenters, pdf1, linewidth=2, label='Simulated PDF')
    plt.plot(binCenters, pdf2, linewidth=2, label='Empirical PDF')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.title(f'{name} with {version}')

    os.makedirs('outputs', exist_ok=True)
    plt.text(binCenters.min() + 0.05*(binCenters.max()-binCenters.min()), max(max(pdf1), max(pdf2))*0.9, f'KS = {ksStat:.3f}', fontsize=12)
    plt.savefig(f"./outputs/{name}_FCD_pdf({version}).png", dpi=300, bbox_inches='tight')
    plt.close()
    return ksStat
