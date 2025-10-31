import numpy as np

#Categories of polymorphic sites taking into account the direction of mutation 

#[in practice, by means of an outgroup and whether the derived state is fixed 

#in population 1 (f1, f1x2) or in population 2 (f2, f2x1)

# Function to count uninterrupted sequences of 0s and 1s

def count_uninterrupted_sequences(sequence):

    count = 1

    current_value = sequence[0]

    modif_sequence = sequence[1:]

    for value in modif_sequence:

        if value != current_value :

            count += 1 

        current_value = value

    return count

def count_distance_sequences (sequence):

    distances = []

    # Initialize the list of groups

    groups = []

    # Initialize variables for the current group

    current_group = []

    current_value = None

    # Iterate over the input list

    for value in sequence:

        # If the current value is different from the previous value,

        # start a new group

        if value != current_value:

            if current_group:

                groups.append(current_group)

                current_group = []

        # Add the current value to the current group

        current_group.append(value)

        # Update the current value

        current_value = value

    # Add the last group to the list of groups

    if current_group:

        groups.append(current_group)


    for group in groups :

        if group[0]== 0:

            distances.append(len(group))

        elif len(group) >= 2 :

            for i in range (len(group)-1):

                distances.append(0)

    return distances

def get_indiv_matrix(haplotype_matrix):

    pop_1 = []

    nbr_indiv_pop = int(haplotype_matrix.shape[1])

    # Sum the columns 2 by 2 for each row

    left_column = haplotype_matrix[:, :nbr_indiv_pop]

    for i in range (0,nbr_indiv_pop,2):

        pop_1.append(np.sum(left_column[:,i:i+2],axis = 1))

    return pop_1

def stat_comp(matrix_comp):

    # Define the set of values to be kept

    values_to_keep = {1, 3, 5, 7}


    # Use numpy's logical_or.reduce to apply multiple conditions

    condition = np.isin(matrix_comp, list(values_to_keep))


    # Initialize an array of zeros with the same shape as matrix_comp

    filtered_array = np.zeros_like(matrix_comp)


    # Set positions to 1 where the condition is true

    filtered_array[condition] = 1

    values_to_change_to_1 = {2, 5, 6, 7, 8}

    Rf = np.isin(matrix_comp, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {4}

    Rs = np.isin(matrix_comp, list(values_to_change_to_1)).astype(int)

    # Define the set of values to be changed to 1

    values_to_change_to_1 = {3, 5}

    # Use vectorized operations with numpy to replace values

    Wx1 = np.isin(matrix_comp, list(values_to_change_to_1)).astype(int)

    Wx1F = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {3}

    Wx1n = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {2,3,5,6,7,8}

    Wx1Fn = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {1,7}

    # Use vectorized operations with numpy to replace values

    Wx2 = np.isin(matrix_comp, list(values_to_change_to_1)).astype(int)

    Wx2F = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {1}

    Wx2n = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    values_to_change_to_1 = {1,2,5,6,7,8} # pas les bons a changer

    Wx2Fn = np.isin(filtered_array, list(values_to_change_to_1)).astype(int)

    return Rs,Rf,Wx1,Wx2,Wx1n,Wx2n,Wx1F,Wx2F,Wx1Fn,Wx2Fn

def W_computation(W_array,seq_len):

    W = 0

    k = seq_len #number of variable sites in the summary sequence

    d = len(W_array) # number of subintervals

    for distance in W_array:

        W+=abs((distance/(k-2)-(1/d))) #distance is the length of the ith sub-interval

    return W/2

def comp_pop_from_haplo (haplo_1, haplo_2):

    ind_haplo_1 = get_indiv_matrix(haplo_1)

    ind_haplo_2 = get_indiv_matrix(haplo_2)

    nbr_indiv_1 = np.shape(ind_haplo_1)[0]

    nbr_indiv_2 = np.shape(ind_haplo_2)[0]

    sum_haplo_1 = np.sum(ind_haplo_1,axis = 0)

    sum_haplo_2 = np.sum(ind_haplo_2,axis = 0)

    fixed_derived_1 = np.isin(sum_haplo_1, nbr_indiv_1*2)

    fixed_ancestral_1 = np.isin(sum_haplo_1, 0)

    mat_to_compare_1 = np.ones_like(fixed_derived_1,dtype = int)

    mat_to_compare_1[fixed_derived_1] = 2

    mat_to_compare_1[fixed_ancestral_1] = 0

    fixed_derived_2 = np.isin(sum_haplo_2, nbr_indiv_2*2)

    fixed_ancestral_2 = np.isin(sum_haplo_2, 0)

    mat_to_compare_2 = np.ones_like(fixed_derived_2,dtype = int)

    mat_to_compare_2[fixed_derived_2] = 2

    mat_to_compare_2[fixed_ancestral_2] = 0

    matrix_comp = 3 * mat_to_compare_1 + mat_to_compare_2

    return matrix_comp

def get_all_dist(matrix_comp):

    # Unpack the output of stat_comp

    Rs, Rf, Wx1, Wx2, Wx1n, Wx2n, Wx1F, Wx2F, Wx1Fn, Wx2Fn = stat_comp(matrix_comp)

    # Define a helper function to compute uninterrupted and distance sequences

    def compute_statistics(arr, count_func):

        return [count_func(sub_arr) for sub_arr in arr]

    # Compute statistics for each required array using the helper function

    Rf_stat = compute_statistics(Rf, count_uninterrupted_sequences)

    Rs_stat = compute_statistics(Rs, count_uninterrupted_sequences)

    Wx2s1_dist = compute_statistics(Wx2, count_distance_sequences)

    Wx1s2_dist = compute_statistics(Wx1, count_distance_sequences)

    Wx1_new_dist = compute_statistics(Wx1n, count_distance_sequences)

    Wx2_new_dist = compute_statistics(Wx2n, count_distance_sequences)

    Wx1F_dist = compute_statistics(Wx1F, count_distance_sequences)

    Wx2F_dist = compute_statistics(Wx2F, count_distance_sequences)

    Wx1F_new_dist = compute_statistics(Wx1Fn, count_distance_sequences)

    Wx2F_new_dist = compute_statistics(Wx2Fn, count_distance_sequences)

    return (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist,Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist)

def get_all_stats(matrix_comp, seq_len):

    # Get all distance-related statistics

    stats = get_all_dist(matrix_comp)


    # Extract the relevant lists from the returned tuple

    (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist,

    Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist) = stats

    # Define a helper function to compute W_computation for a list of distances

    def compute_W_statistics(distances, seq_len):

        return [W_computation(dist, seq_len) for dist in distances]

    # Compute W statistics for each distance list

    Wx2s1_stat = compute_W_statistics(Wx2s1_dist, seq_len)

    Wx1s2_stat = compute_W_statistics(Wx1s2_dist, seq_len)

    Wx1F_stat = compute_W_statistics(Wx1F_dist, seq_len)

    Wx2F_stat = compute_W_statistics(Wx2F_dist, seq_len)

    Wx1_new_stat = compute_W_statistics(Wx1_new_dist, seq_len)

    Wx2_new_stat = compute_W_statistics(Wx2_new_dist, seq_len)

    Wx1F_new_stat = compute_W_statistics(Wx1F_new_dist, seq_len)

    Wx2F_new_stat = compute_W_statistics(Wx2F_new_dist, seq_len)

    return (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist, Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist, 
    Wx2s1_stat, Wx1s2_stat, Wx1F_stat, Wx2F_stat, Wx1_new_stat, Wx2_new_stat, Wx1F_new_stat, Wx2F_new_stat)

def get_stat_from_haplo (haplo_1, haplo_2):

    new_haplo_1 = get_indiv_matrix(haplo_1)

    nbr_indiv = len(new_haplo_1)

    new_haplo_2 = get_indiv_matrix(haplo_2)

    if nbr_indiv > len (new_haplo_2):

        nbr_indiv = len (new_haplo_2)

    new_haplo_1 = np.array(new_haplo_1 [:nbr_indiv])

    new_haplo_2 = np.array(new_haplo_2 [:nbr_indiv])

    matrix_comp = 3 * new_haplo_1 + new_haplo_2

    seq_len = matrix_comp.shape[1]

    (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist, Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist, 
            Wx2s1_stat, Wx1s2_stat, Wx1F_stat, Wx2F_stat, 
            Wx1_new_stat, Wx2_new_stat, Wx1F_new_stat, Wx2F_new_stat) = get_all_stats(matrix_comp, seq_len)

    return (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist, 
    Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist, 
    Wx2s1_stat, Wx1s2_stat, Wx1F_stat, Wx2F_stat, 
    Wx1_new_stat, Wx2_new_stat, Wx1F_new_stat, Wx2F_new_stat)

def get_stat_pop_from_haplo (haplo_1, haplo_2):

    matrix_comp = np.array([comp_pop_from_haplo(haplo_1, haplo_2)])

    seq_len = matrix_comp.shape[1]

    (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist, 
    Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist, 
    Wx2s1_stat, Wx1s2_stat, Wx1F_stat, Wx2F_stat, 
    Wx1_new_stat, Wx2_new_stat, Wx1F_new_stat, Wx2F_new_stat) = get_all_stats(matrix_comp, seq_len)

    return (Rf_stat, Rs_stat, Wx2s1_dist, Wx1s2_dist, Wx1F_dist, Wx2F_dist, 
    Wx1_new_dist, Wx2_new_dist, Wx1F_new_dist, Wx2F_new_dist, 
    Wx2s1_stat, Wx1s2_stat, Wx1F_stat, Wx2F_stat, 
    Wx1_new_stat, Wx2_new_stat, Wx1F_new_stat, Wx2F_new_stat)
