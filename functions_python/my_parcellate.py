import numpy as np

def my_parcellate(data, hemisphere):
    parc_name = 'Glasser360'
    parc_file = f'/home/zwang/storage/my/BrainEigenmodes_legacy/data/parcellations/fsLR_32k_{parc_name}-{hemisphere}.txt'
    parc = np.loadtxt(parc_file)
    data_parc = calc_parcellate(parc, data)
    return np.nan_to_num(data_parc, nan=0.0)


def calc_parcellate(parc, data_input):
    num_vertices = parc.shape[0]
    print(num_vertices)
    parcels = np.unique(parc[parc > 0])
    num_parcels = len(parcels)

    if data_input.shape[0] != num_vertices:
        data_input = data_input.T

    data_parcellated = np.zeros((num_parcels, data_input.shape[1]))
    for i, parcel_interest in enumerate(parcels):
        ind_parcel = np.where(parc == parcel_interest)[0]
        data_parcellated[i, :] = np.mean(data_input[ind_parcel, :], axis=0)
    return data_parcellated
