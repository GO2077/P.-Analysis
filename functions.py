import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors



# Determine matches

def cat(ped, snp):
    # Check the PEDSEX and SNPSEX values against one another and return a value for the possible conditions
    # 0=PROBLEM; 1=PlinkError; 2=SexCheckCorrect; 3=SexCheckIncorrect


    if ped == 0:
        return 0
    elif snp == 0:
        return 1
    elif ped == snp:
        return 2
    else:
        return 3



# Rotate

def rotate(mat, fac):
    # Take the rows and columns of the input matrix and compute a 90° rotation

    case_layer = mat.xs(fac, axis=1, level='Field')
    pos_layer  = mat.xs('code_position', axis=1, level='Field')

    rows, cols = case_layer.shape

    case_flipped = np.empty((cols, rows), dtype=object)
    pos_flipped = np.empty((cols, rows), dtype=object)

    for i in range(rows):
        for j in range(cols):
            case_flipped[j, rows - 1 - i] = case_layer.iloc[i, j]
            pos_flipped[j, rows - 1 - i] = pos_layer.iloc[i, j]

    new_columns = pd.MultiIndex.from_product([case_layer.index, [fac, 'code_position']],
                                             names=mat.columns.names)

    combined_values = np.empty((case_flipped.shape[0], case_flipped.shape[1]*2), dtype=object)
    combined_values[:, ::2] = case_flipped
    combined_values[:, 1::2] = pos_flipped

    new_index = case_layer.columns

    flipped_df = pd.DataFrame(combined_values, index=new_index, columns=new_columns)
    return flipped_df



# Flip

def flip(mat, fac):
    # Take the rows and columns of the input matrix and transpose along the diagonal

    case_layer = mat.xs(fac, axis=1, level='Field')
    pos_layer  = mat.xs('code_position', axis=1, level='Field')

    rows, cols = case_layer.shape

    case_flipped = np.empty((cols, rows), dtype=object)
    pos_flipped = np.empty((cols, rows), dtype=object)

    for i in range(rows):
        for j in range(cols):
            case_flipped[j, i] = case_layer.iloc[i, j]
            pos_flipped[j, i] = pos_layer.iloc[i, j]

    new_columns = pd.MultiIndex.from_product([case_layer.index, [fac, 'code_position']],
                                             names=mat.columns.names)

    combined_values = np.empty((case_flipped.shape[0], case_flipped.shape[1]*2), dtype=object)
    combined_values[:, ::2] = case_flipped
    combined_values[:, 1::2] = pos_flipped

    new_index = case_layer.columns

    flipped_df = pd.DataFrame(combined_values, index=new_index, columns=new_columns)
    return flipped_df




# Shift

def shift(mat, fac):
    # Reshift the plate while retaining the original dimensions

    case_layer = mat.xs(fac, axis=1, level='Field')
    pos_layer  = mat.xs('code_position', axis=1, level='Field')

    rows, cols = case_layer.shape

    shift_case = case_layer.values.flatten(order='C')
    shift_pos= pos_layer.values.flatten(order='C')

    if len(shift_case) % rows != 0:
        raise ValueError(f"Total elements ({len(shift_case)}) is not divisible by amount of rows ({rows})."
                         " Recheck input.")
    n_cols = len(shift_case) // rows

    reshaped_case = shift_case.reshape((rows, cols), order='F')
    reshaped_pos = shift_pos.reshape((rows, cols), order='F')

    new_columns = pd.MultiIndex.from_product([range(1, cols+1), [fac, 'code_position']],
                                             names=mat.columns.names)


    combined_values = np.empty((rows, cols*2), dtype=object)
    combined_values[:, ::2] = reshaped_case
    combined_values[:, 1::2] = reshaped_pos

    row_labels = [chr(ord('A') + i) for i in range(rows)]

    shift_matrix = pd.DataFrame(combined_values, index=row_labels, columns=new_columns)
    return shift_matrix



# Visualize

def visualize_matrix(matrix, title, fac):
    matrix_layer = matrix.xs(fac, axis=1, level='Field')
    matrix_index = matrix.loc[:, pd.IndexSlice[:, 'code_position']].values

    matrix_numeric = pd.to_numeric(matrix_layer.stack(), errors='coerce').unstack().fillna(-1).values

    cmap = mcolors.ListedColormap(['gray', 'green', 'blue', 'orangered'])
    #bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    X, Y = np.meshgrid(np.arange(matrix_numeric.shape[1]+1), np.arange(matrix_numeric.shape[0]+1))

    plt.figure(figsize=(6,6))
    plt.pcolormesh(X, Y, matrix_numeric, cmap=cmap, norm=norm, edgecolors='black', linewidth=1)

    plt.xticks(np.arange(0.5, matrix_numeric.shape[1], 1), labels=matrix_layer.columns.get_level_values(0))
    plt.yticks(np.arange(0.5, matrix_numeric.shape[0], 1), labels=matrix_layer.index)

    plt.gca().invert_yaxis()

    for i in range(matrix_numeric.shape[0]):
        for j in range(matrix_numeric.shape[1]):
            val = matrix_numeric[i, j]
            sign = matrix_index[i, j]
            display_sign = '' if sign is None else str(sign)
            plt.text(j + 0.5, i + 0.5, display_sign, ha='center', va='center', color='black')

    plt.title(title)
    plt.show()
