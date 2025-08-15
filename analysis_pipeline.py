import numpy as np
import matplotlib.pyplot as plt
import os

# Robust relative imports to allow running from package or as script
try:
    from .functions_python.draw_tc_samples import draw_tc_samples
    from .functions_python.my_parcellate import my_parcellate
    from .functions_python.draw_FC import draw_FC
    from .functions_python.draw_FC_corr_scatter import draw_FC_corr_scatter
    from .functions_python.draw_FCD_sliding import draw_FCD_sliding
    from .functions_python.draw_FCD_phase import draw_FCD_phase
    from .functions_python.draw_FCD_pdf_KS import draw_FCD_pdf_KS
except Exception:
    from functions_python.draw_tc_samples import draw_tc_samples
    from functions_python.my_parcellate import my_parcellate
    from functions_python.draw_FC import draw_FC
    from functions_python.draw_FC_corr_scatter import draw_FC_corr_scatter
    from functions_python.draw_FCD_sliding import draw_FCD_sliding
    from functions_python.draw_FCD_phase import draw_FCD_phase
    from functions_python.draw_FCD_pdf_KS import draw_FCD_pdf_KS


def analysis_pipeline(bold_full, name, T, cortex, surface=None):
    cortex_ind = np.where(cortex)[0]

    os.makedirs('./outputs', exist_ok=True)

    mean_activity = np.mean(bold_full, axis=0)
    plt.figure(figsize=(10, 6))
    plt.plot(mean_activity, linewidth=1.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Mean BOLD')
    plt.title('Mean BOLD Across Vertices')
    plt.grid(True)
    plt.savefig(f'./outputs/{name}_mean_BOLD.png', dpi=300, bbox_inches='tight')
    plt.close()

    draw_tc_samples(bold_full, cortex_ind, f"{name}_BOLD_TC_samples", T, 8, False)

    data_parc = my_parcellate(bold_full, 'lh')
    print("Size of parcelled bold:")
    print(data_parc.shape)

    FC_vec = draw_FC(data_parc, name)

    FC_corr = {
        'Minimal': draw_FC_corr_scatter(FC_vec, name, 'Minimal'),
        'ICA': draw_FC_corr_scatter(FC_vec, name, 'ICA'),
        'ICA_GSR': draw_FC_corr_scatter(FC_vec, name, 'ICA+GSR'),
        'ICA_GSR_WM_CSF_MT_CEN': draw_FC_corr_scatter(FC_vec, name, 'ICA+GSR+WM_CSF_MT_CEN'),
    }

    FCD_sliding = draw_FCD_sliding(data_parc, 83, name)
    KS = {
        'sliding_ICA': draw_FCD_pdf_KS(FCD_sliding, "sliding_ICA", name),
    }

    FCD_phase = draw_FCD_phase(data_parc, name)
    KS.update({
        'phase_Minimal': draw_FCD_pdf_KS(FCD_phase, "phase_Minimal", name),
        'phase_ICA': draw_FCD_pdf_KS(FCD_phase, "phase_ICA", name),
    })

    return {'FC_corr': FC_corr, 'KS': KS}

