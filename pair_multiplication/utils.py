# utils.py
from collections import Counter
import numpy as np
import itertools as it

def _Nc_positive(Nc:int):
        
    Nc = max([Nc,0])
        
    return Nc
    
def partition_tuplify(part):

    if not isinstance(part, tuple):
        part = (part,)
        
    return part

def is_monotonic_decreasing(lis):
    if len(lis)>0:
        return all(x >= y for x, y in zip(lis, lis[1:]))
    else:
        return True
def is_monotonic_increasing(lis):
    if len(lis)>0:
        return all(x <= y for x, y in zip(lis, lis[1:]))
    else:
        return True
def make_tuple(lis):
    
    return tuple(Counter(lis).values())


def is_serial_increasing(lis):
    if len(lis)>0:
        return list(set(lis)) == list(np.arange(1,np.max(lis)+1))
    else:
        return True
def is_serial_increasing_from_its_min(lis):
    if len(lis)>0:
        return list(set(lis)) == list(np.arange(np.min(lis),np.max(lis)+1))
    else:
            return True
def is_serial_decreasing(lis):
    if len(lis)>0:
        return list(set(lis)) == list(np.flip(np.arange(1,np.max(lis)+1)))
    else:
        return True

def is_serial_decreasing_from_its_min(lis):
    
    if len(lis)>0:
        return list(set(lis)) == list(np.flip(np.arange(np.min(lis),np.max(lis)+1)))
    else:
        return True


def is_subword(lis):
    
    counter_lis = list(Counter(lis).values())
    counter_keys = list(Counter(lis).keys())
    return is_monotonic_decreasing(counter_lis) and \
        (sorted(counter_keys) == list(range(1, int(max(counter_keys))+1)))

def is_word(lis):
    
    return np.all([is_subword(lis[0:indx+1]) for indx in range(len(lis)-1)])

def is_column_increasing(lists):
    
    return np.array([np.all([is_serial_increasing(li) for li in 
                        list({re:lis.imag[lis.real.astype(int)==re.astype(int)] for re in lis.real}.values())])
                     for lis in lists])

def is_column_increasing_from_its_min(lists):

        return np.array([np.all([is_serial_increasing_from_its_min(li) for li in 
                        list({re.astype(int):lis.imag.astype(int)[lis.real.astype(int)==re.astype(int)] 
                              for re in lis.real}.values())])
                     for lis in lists])

def is_column_decreasing_from_its_min(lists):

        return np.array([np.all([is_serial_decreasing_from_its_min(li) for li in 
                        list({re:lis.imag[lis.real.astype(int)==re.astype(int)] for re in lis.real}.values())])
                     for lis in lists])

def is_column_decreasing(lists):

        return np.array([np.all([is_serial_decreasing(li) for li in 
                        list({re:lis.imag[lis.real.astype(int)==re.astype(int)] for re in lis.real}.values())])
                     for lis in lists])

def is_word_column_admissible(lis):
    
    if len(lis)>0:
        imag_part = np.all([is_serial_increasing(li) for li in 
                            list({re:lis.imag[lis.real==re] for re in lis.real}.values())])
        real_counter = list(Counter(lis.real).values())
        real_entries = list(Counter(lis.real).keys())
        real_part = is_monotonic_decreasing(real_counter) and\
                    np.all(set(np.unique(real_entries))==set(np.arange(1,max(real_entries)+1)))

        return real_part and imag_part 
    else:
        return True


def get_all_combinations(word, num_letters):
    
    all_combs = np.array(list(it.combinations(word, num_letters)))

    return all_combs

def get_all_permutations(word,num_letters):
    
    all_perms = np.array(list(it.permutations(word, num_letters)))
    
    return all_perms

def remove_all_entries(word, combs):
    
    word_counted = Counter(word)
    
    combs_popped = np.array([tuple((word_counted-Counter(list(combs[ind]))).keys())
                      for ind in range(len(combs))])
    
    return combs_popped

def check_row_admissible(removingB, removingA):
    
    row_counter = Counter(removingB.real)
    
    return np.all([len(dict(zip(removingA.imag[removingB.real==row], 
                         removingA.real[removingB.real==row]))) == row_counter[row]
                   for row in list(row_counter.keys())])

def check_column_admissible(cols, entries, num_letters):
    
    return len(list(set([(cols[indi],entries[indi]) for indi in range(num_letters)]))) >= num_letters


def sort_word_index_LR(complex_list):
    
    return sorted(range(len(complex_list)), key=lambda i: (complex_list[i].real, -complex_list[i].imag))

def sort_word(complex_list):
    
    complex_list = np.array(complex_list)
    
    return np.array(sorted(complex_list, key=lambda x: (x.real, x.imag)))


def check_monotonic_increasing(row,mask,real_part):
    filtered_row = np.trim_zeros(real_part[row][mask[row]],trim='b')
    return np.all(np.diff(filtered_row) >= 0) and not 0 in filtered_row
def check_monotonic_increasing_strict(col,mask,real_part):
    filtered_row = np.trim_zeros(real_part[:,col][mask[:,col]],trim='b')
    return np.all(np.diff(filtered_row) > 0) and not 0 in filtered_row

def is_monotonic_increasing(arr):
    # Get the real and imaginary parts of the array
    real_part = arr.real
    imag_part = arr.imag
    
    # Create a mask for entries with non-zero imaginary part
    mask = imag_part == 0  # True where imaginary part is zero, False otherwise
    #print(mask)
    # Function to check monotonicity in a 1D array with a mas
    # Check each row and column for monotonic increasing real parts
    rows_increasing = all(check_monotonic_increasing(row,mask,real_part) for row in range(arr.shape[0]))
    cols_increasing = all(check_monotonic_increasing_strict(col,mask,real_part) for col in range(arr.shape[1]))

    # Result is True only if both rows and columns are monotonic increasing
    return rows_increasing and cols_increasing


def make_2d_array(perm):
    
    arr = np.ones((len(perm),perm[0]))
    
    for row in range(len(perm)):
        
        arr[row][perm[row]:] = 0*arr[row][perm[row]:]
        arr[row][0:perm[row]] = (row+1)*arr[row][0:perm[row]]
        
    return arr


def LR_multiple(permA,permB):
    
    
    diagB = make_2d_array(permB)
    map_diag = np.ones((len(permA)+len(permB),permA[0]+permB[0]),dtype=complex)*0
    map_diag[0:len(permB),0:permB[0]] = diagB*1j
    
    which_rows_to_add = [np.arange(row+1,len(permB)+row+2) for row in range(len(permA))]

    locs = [np.array(list(it.combinations_with_replacement(which_rows_to_add[row],permA[row]))) 
        for row in range(len(permA))]
    
    saved_cands = [np.array(map_diag)]

    for rowa in range(len(permA)):

        intermittent_cands = []

        for scand in saved_cands:
            for comb in locs[rowa]:    

                cand = np.array(scand)
                for row_lab in list(set(comb)):

                    row = comb[np.where(comb==row_lab)[0]]
                    r=row_lab-1
                    start=np.where(cand[r]==0)[0][0]
                    cand[r,start:start+len(row)] = np.ones(len(row))*(rowa+1)

                if is_monotonic_increasing(cand):
                    intermittent_cands.append(cand)
        saved_cands = intermittent_cands

    saved_cands = np.array(saved_cands)
    
    flipped_cands = np.flip(saved_cands.real,axis = -1)
    flattened = np.array([arr.flatten() for arr in flipped_cands])
    keep = [is_word(arr[arr>0]) and len(arr[arr>0])==sum(permA) for arr in flattened]
    kept_cands = saved_cands[keep]
    
    results = [tuple([np.max(np.where(cnd!=0)[0])+1 for cnd in cand if np.any(cnd!=0)]) for cand in kept_cands]
    
    return dict(Counter(results))
    

def superpose(num_letters, wordB,wordA):

    all_supA = get_all_permutations(wordA,num_letters)
    all_supB = get_all_combinations(wordB,num_letters)

    all_remA = remove_all_entries(wordA,all_supA)
    all_remB = remove_all_entries(wordB,all_supB)

    keepB = [is_word_column_admissible(all_remB[ind]) 
             for ind in range(len(all_remB))]
    keepA = np.all([[is_word_column_admissible(all_remA[ind])
                     for ind in range(len(all_remA))],
                    is_column_increasing(all_remA),
                   ], axis =0)
    
    
    candsB = all_remB[keepB]
    supB = all_supB[keepB]

    candsA = all_remA[keepA]
    supA = all_supA[keepA]

    uniqueA = np.array([np.where(np.sum(abs(supA.real-el), axis = -1)==0)[0][0] 
                        for el in list(np.unique(supA.real,axis = 0))])

    candsA_unq = candsA[uniqueA]
    supA_unq = supA[uniqueA]
    
    multiplicity_counter = [
                     str([sorted(supA_unq.real[inda][supB[indb].real==len(list(set(wordB.real)))-row])
                      for row in range(len(list(set(wordB.real))))])
                    for indb in range(len(candsB)) for inda in range(len(candsA_unq))
                    if is_word(np.concatenate((candsA_unq[inda].real,
                                               supA_unq[inda].real,
                                              )))
                    and check_column_admissible(supB[indb].imag,
                                               supA_unq[inda].real,
                                               num_letters)
                    and check_row_admissible(supB[indb],
                                             supA_unq[inda]
                                            )
                    ]
    
    mults_cands = [(make_tuple(candsB[indb].real),make_tuple(candsA_unq[inda].real))
                for indb in range(len(candsB)) for inda in range(len(candsA_unq))
                if is_word(np.concatenate((candsA_unq[inda].real,
                                           supA_unq[inda].real,
                                          )))
                and check_column_admissible(supB[indb].imag,
                                           supA_unq[inda].real,
                                           num_letters)
                and check_row_admissible(supB[indb],
                                         supA_unq[inda]
                                        )
                ]
    
    mults = dict(Counter(list(dict(zip(multiplicity_counter,mults_cands)).values())))

    return mults
    
    
            
def compose_pair(pair):
    
    tmax = min([sum(pair.pair[0].n),sum(pair.pair[1].n)])
    
    wordA = pair.pair[1].word
    wordB = 0
    
    all_combinations = list(it.chain(*[superpose(t,pair) for t in range(tmax)]))
    
    all_combinations_admissible = drop_invalid(all_combinations)
    
    
def extend_partition(perm,new_length):
    
    perm = list(perm)
    perm += list(np.zeros(new_length-len(perm)))
    
    return perm
