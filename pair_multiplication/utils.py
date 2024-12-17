# utils.py
import warnings
from collections import Counter
import numpy as np
import itertools as it

def _Nc_positive(Nc:int):
        
    Nc = max([Nc,0])
        
    return Nc
    
def partition_tuplify(part):

    if not isinstance(part, tuple):
        if 0==part:
            return ()
    
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
    
    if len(lis)>1:
        checker = [is_subword(lis[0:indx+1]) for indx in range(len(lis))]
        return np.all(checker)
    elif len(lis)==1:
        if max(lis)<=1:
            return True
        else:
            return False
    else:
        return True

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

def is_monotonic_increasing_inarray(arr):
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

def make_2d_array_ones(perm):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        arr = np.zeros((1,1))
        if len(perm)>0:
            arr = np.ones((len(perm),perm[0]))
            
            for row in range(len(perm)):
                
                arr[row][perm[row]:] = 0*arr[row][perm[row]:]
            
        return arr


def make_2d_array_LR(perm):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        arr = np.zeros((1,1))
        if len(perm)>0:
            arr = np.ones((len(perm),perm[0]))
            
            for row in range(len(perm)):
                
                arr[row][perm[row]:] = 0*arr[row][perm[row]:]
                arr[row][0:perm[row]] = (row+1)*arr[row][0:perm[row]]
            
        return arr


def LR_multiple(permA,permB):

    if len(permA)==0:
        return {permB:1},[permB],np.array([[0]])
    elif len(permB)==0:
        return {permA:1},[permA],np.array([[0]])
    else:
        diagB = make_2d_array_LR(permB)
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

                    if is_monotonic_increasing_inarray(cand):
                        intermittent_cands.append(cand)
            saved_cands = intermittent_cands

        saved_cands = np.array(saved_cands)
        
        flipped_cands = np.flip(saved_cands.real,axis = -1)
        flattened = np.array([arr.flatten() for arr in flipped_cands])
        keep = [is_word(arr[arr>0]) and len(arr[arr>0])==sum(permA) for arr in flattened]
        kept_cands = saved_cands[keep]
        
        results = [tuple([np.max(np.where(cnd!=0)[0])+1 for cnd in cand if np.any(cnd!=0)]) for cand in kept_cands]
        
        return dict(Counter(results)), results, kept_cands
    

def superpose(num_letters, wordB,wordA):

    if np.all(wordA==0j) and np.all(wordB ==0j):
        perm = ((),())
        return {perm:1},[perm],['[]']
    elif np.all(wordA==0j):
        perm = tuple([sum(wordB.real==row+1) for row in range(max(wordB.real).astype(int))])
        return {(perm,()):1},[(perm,())],['[]']
    elif np.all(wordB==0j):
        perm = tuple([sum(wordA.real==row+1) for row in range(max(wordA.real).astype(int))])
        return {((),perm):1},[((),perm)],['[]']
    else:
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

        return mults, mults_cands, multiplicity_counter
    

    
def extend_partition(perm,new_length):
    
    perm = list(perm)
    perm += list(np.zeros(new_length-len(perm)))
    
    return perm
    
    
##################################### For Column Adding Methods #########################################
def col_add(brA,brB):
    col = 0
    if len(brA)>0:
        col += brA[0]
    if len(brB)>0:
        col += brB[0]
        
    return col

def make_2d_array(perm, column_sort=False, dims = None):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        dim = ((1,1))
        if not dims is None:
            dim = dims
        
        arr = np.zeros(dim, dtype=np.complex128)
        if column_sort:
            arr = arr*(0.5+0.5j)
        if len(perm)>0 or np.any(dim)>1:
            
            rows = np.max([dim[0],len(perm)])
            cols = np.max([dim[1],perm[0]])
            
            arr = np.zeros((rows,cols),dtype=complex)
            
            for row in range(len(perm)):
                
                arr[row][perm[row]:] = 0*arr[row][perm[row]:]
                if column_sort:
                    arr[row][0:perm[row]] = np.arange(1,perm[row]+1)#(row+1)*arr[row][0:perm[row]]
                else:
                    arr[row][0:perm[row]] = np.ones(perm[row])*(0.5+0.5j)
        
        return arr
    
def make_pair_array(pair,column_sort = False):
    
    part_br = pair[0]
    part_ubr = pair[1]

    arr_ubr = make_2d_array(part_ubr,column_sort=column_sort)
    arr_br = make_2d_array(part_br,column_sort=column_sort) 
    
    return arr_br*-1j,arr_ubr
    
def make_candidate_array(pairB,pairA):
    
    brB = pairB[0]
    brA = pairA[0]
    
    ubrB = pairB[1]
    ubrA = pairA[1]
    
    rows_br = len(brB)+len(brA)
    cols_br = col_add(brA,brB)
    
    rows_ubr = len(ubrB)+len(ubrA)
    cols_ubr = col_add(ubrA,ubrB)
    
    cols = np.max([cols_br,cols_ubr])
    
    arr_ubr = make_2d_array(ubrB,dims = (rows_ubr+rows_br,cols), column_sort=False)
    arr_br = make_2d_array(brB,dims = (rows_br+rows_ubr,cols), column_sort=False)*-1 
    
    return arr_br, arr_ubr

def get_entries_to_add(array,t):
    
    arr = array.copy()#np.array(array, dtype=np.complex128)
    coast_arr = arr.real==0
    coast = np.argwhere(coast_arr==True)

    if len(coast) > 0:
    
        ents = np.array([cl
                         for cl in coast 
                     if (cl[0] < np.min((coast[:,0])[coast[:,1]==cl[1]])+t)])
    
        return ents
    else:
        return np.array([])


def get_entries_to_remove(array,t):
    
    arr = array.copy()#np.array(array)
    coast_arr = (arr.imag.astype(np.int64)==0) * (arr.real!=0)
    coast = np.argwhere(coast_arr==True)
    if len(coast)>0:
        ents = np.array([cl
                         for cl in coast 
                     if (cl[0] >= np.max((coast[:,0])[coast[:,1]==cl[1]])-t+1)])
        return ents
    else:
        return np.array([])

def in_different_rows(cands):
    try:
        cands_kept = np.array([cand for cand in cands if len(list(set(cand[:,0])))==len(cands[0,:,0])
                      ])
    
        return cands_kept
    except Exception as e:
        return np.array([[]])
        
def check_words(u_br,u_nbr,min_n0,return_n0=False):

    concd = get_overlap(u_br,u_nbr,min_n0)
    
    rtn = False
    if not (concd is None):
        rtn = True
    if return_n0:
        
        if rtn:
            br_word_array = np.flip(concd,axis = 0)
            ubr_word_array = np.flip(concd.T,axis = 1)

            n0 = len(ubr_word_array[0])
            min1s = np.sum(ubr_word_array.real == -1,axis = 0)+\
                    np.sum(ubr_word_array.imag == -1,axis = 0)

            ones = np.sum(ubr_word_array.real == 1,axis = 0)+\
                   np.sum(ubr_word_array.imag == 1,axis = 0)

            step_up = [0]+[np.sum(ones[0:k+1])-np.sum(min1s[0:k+1]==0) for k in range(len(ones)) 
                          ]

            n0 += max(step_up)

            n0_cont = np.zeros((1,len(ubr_word_array[0])))
            n0_cont[-1] = n0
            checker = np.concatenate((ubr_word_array,n0_cont),axis = 0)

            return (n0,checker,rtn)
        else:
            return (0,np.array([[]], dtype=np.int64),False)
        
    return rtn


def check_shape(diag):
    
    mask = (diag.real!=0) * (diag.imag.astype(np.int64)<=0)

    dig = np.zeros(diag.shape, dtype=np.int64)
    dig[mask] += 1
    
    perm = np.trim_zeros(np.sum(dig,axis = -1).astype(np.int64))
    
    gaps = np.all(np.apply_along_axis(is_monotonic_decreasing,0,dig))
    
    return is_monotonic_decreasing(perm)*gaps

def check_is_diagram(diag, num_entries, val):
    
    mask = (diag.real!=0) * (diag.imag.astype(np.int64)<=0)

    dig = np.zeros(diag.shape, dtype=np.int64)
    dig[mask] += 1
    
    perm = np.trim_zeros(np.sum(dig,axis = -1).astype(np.int64))
    
    check_re = np.sum(diag.real == val.real) == num_entries
    check_im = np.sum(diag.imag == val.imag) == num_entries
    
    checker_num_added = check_re+check_im
    
    return check_shape(diag) * check_shape(diag.T) * checker_num_added


def replace_column(diag, adding_col, label):
    
    diag = diag.copy()
    rows,cols = adding_col.T
    diag[rows,cols] = label
    return diag


def add_column(diag, adding_col, label):
    
    diag = diag.copy()
    rows,cols = adding_col.T
    
    diag[rows,cols] = diag[rows,cols].real.astype(np.int64)+1j*diag[rows,cols].imag.astype(np.int64)
    
    diag[rows,cols] += label
    
    return diag


def columns_align(cbr,cubr):

    c0_cond = np.all((cbr.real!=cubr.imag)[cubr.imag.astype(np.int64)<0])*\
              np.all((cbr.imag!=cubr.real)[cubr.real.astype(np.int64)>0])*\
              np.all((abs(cbr.real)!=0.5)*(abs(cubr.real)!=0.5))
    
    if c0_cond:
        c1_cond = ((cubr*1j - cbr).imag)[(cubr.real > 0)*(cbr.imag > 0)]

        c2_cond = (cbr.real - cubr.imag)[(cubr.imag < 0)*(cbr.real < 0)]
        
        c3_cond = (~((cubr.real != 0)*(cbr.imag == 0)))[(cubr.imag < 0)*(cbr.real < 0)]
        
        align = False
        
        if len(c1_cond)==len(cbr):
        
            align = (np.all(c1_cond.astype(np.int64)>0))*(len(c1_cond)==len(cbr))*np.all((c2_cond.astype(np.int64)>0))
            
        elif len(c2_cond)==len(cbr) and np.all(c1_cond>0):
            align = (np.all(c1_cond.astype(np.int64)>0))*(len(c1_cond)==len(cbr))*np.all((c2_cond.astype(np.int64)>0))+\
                (np.all((c2_cond.astype(np.int64)>0))*(len(c2_cond)==len(cbr))*np.all(c3_cond))
        
        
        return align
    else:

        return False
        
        
def get_full_diag(br,ubr,arg_min):
    
    concd = np.concatenate((np.flip(np.flip(br[0:arg_min],axis = 1),axis =0),ubr[0:arg_min]), axis = 1)
    
    return concd
    
    
def check_wording(concd):
    
    br_word_array = np.flip(concd,axis = 0)
    mask_br = (br_word_array.real.astype(np.int64) < 0) + (br_word_array.imag.astype(np.int64) < 0)
    br_w = br_word_array[mask_br]
    
    br_word = np.zeros(len(br_w))

    br_word[br_w.real<0] -= br_w.real[br_w.real<0]
    br_word[br_w.imag<0] -= br_w.imag[br_w.imag<0]
    
    
    ubr_word_array = np.flip(concd.T,axis = 1)
    
    mask_ubr = (ubr_word_array.imag.astype(np.int64) > 0) + (ubr_word_array.real.astype(np.int64) > 0)
    ubr_w = ubr_word_array[mask_ubr]

    ubr_word = np.zeros(len(ubr_w))

    ubr_word[ubr_w.imag>0] += ubr_w.imag[ubr_w.imag>0]
    ubr_word[ubr_w.real>0] += ubr_w.real[ubr_w.real>0]
    
    rtn = is_word(br_word) and is_word(ubr_word)
    return rtn
    
def get_overlap(br,ubr,min_n0):
    
    first_cols=np.array([br[:,0], ubr[:,0]])
    arg_min = np.sum(np.argmin(abs(first_cols),axis = 1),axis = 0)
    
    concd = np.concatenate((np.flip(np.flip(br[0:arg_min],axis = 1),axis =0),ubr[0:arg_min]), axis = 1)
    
    c_ubr = np.flip(np.trim_zeros(first_cols[1],'b'))

    c_br = np.flip(np.trim_zeros(first_cols[0],'b'))

    kmax = min(len(c_ubr),len(c_br))

    overlap = max([0]+[k+1 for k in range(kmax)
            if columns_align(np.flip(c_br[0:k+1]),c_ubr[0:k+1]) and
            check_wording(get_full_diag(br,ubr,arg_min-(k+1)))
                      ])
    
    arg_min -= overlap
    concd = get_full_diag(br,ubr,max(arg_min,min_n0))
    
    if not check_wording(concd):
        concd = None
    
    return concd
    


def remove_combos_with_insufficient_gaps(removal_combos,diag):
    
    gaps_br = np.sum((diag.real == -1),axis =-1)
    
    recombs = np.array([removal_combos[comb] for comb in range(len(removal_combos))
        if np.all([(np.sum(gaps_br[0:removal_combos[comb,ind,0]+1]==0) >= len(removal_combos[comb,0:ind+1,0]))
                   for ind in range(len(removal_combos[comb]))]) 
                       ], dtype=np.int64)
    
    if len(recombs)==0:
        return np.array([[]], dtype=np.int64)
    
    return recombs
    
    
def count_steps_up(cad,crm, nc_br,nc_ubr):
    
    zer_ind = np.where(np.sum(nc_br!=0,axis = -1)==0)[0]
    zind = len(nc_br)
    if len(zer_ind)>0:
        zind = zer_ind[0]
    gaps_ubr = np.sum((nc_ubr.imag == -1)*(nc_ubr != 0),axis =-1)[nc_ubr[:,0]!=0]
    gaps_br = np.sum(np.sum((nc_br.real == -1)*(nc_br.imag!=1),axis =-1)[0:zind]==0)
    
    steps = np.array([1 for ind in range(len(cad))
                      if np.sum(gaps_ubr[cad[ind,0]:]==0) < ind+1], dtype=np.int64)
    
    return np.max([np.sum(steps)-gaps_br,0])
    
def unique_permutations_with_N0s(diags_br,diags_ubr,n0s):

    min_n0 = n0s
    
    mask_ubr = np.array(diags_ubr.real)!=0
    ubrs = np.zeros(diags_ubr.shape)
    ubrs[mask_ubr] += 1
    
    perms_ubr = np.sum(ubrs,axis = -1)
    
    mask_br = (diags_br.real!=0) * (diags_br.imag.astype(np.int64)==0)
    brs = np.zeros(diags_br.shape)
    brs[mask_br] += 1
    
    perms_br = np.sum(brs,axis = -1)

    pairs = [(tuple(np.trim_zeros(perms_br[ind]).astype(np.int64)),
              tuple(np.trim_zeros(perms_ubr[ind]).astype(np.int64)))
             for ind in range(len(perms_br))]
    
    return pairs, n0s.astype(int)
    
def pair_multipy(pairB,pairA):
    
    brA,ubrA = make_pair_array(pairA,column_sort=True)
    brB,ubrB = make_pair_array(pairB)

    diagC_br, diagC_ubr = make_candidate_array(pairB,pairA)
    
    min_n0 = max(len(pairA[1])+len(pairA[0]),len(pairB[0])+len(pairB[1]))
    
    checking_diags = []
    candidates_br = [diagC_br]

    candidates_ubr = [diagC_ubr]

    cols_br = np.array(list(set(-brA.imag.astype(int).flatten())))
    cols_br = np.flip(cols_br[cols_br>0])
    
    for col in cols_br:

        num_boxes = np.sum(-brA.imag.astype(int)==col)

        cands_br = []
        cands_ubr = []

        for cind in range(len(candidates_br)):
            c_br = candidates_br[cind]
            c_ubr = candidates_ubr[cind]
            tmax = min([np.sum(c_ubr.real!=0),num_boxes])
            
            for t in range(tmax+1):

                nmint = num_boxes - t
                where_to_add = get_entries_to_add(c_br,nmint)
                where_to_remove = get_entries_to_remove(c_ubr,t)

                addition_combos = in_different_rows(np.array(list(it.combinations(where_to_add, nmint))))
                removal_combos = in_different_rows(np.array(list(it.combinations(where_to_remove, t))))

                for cad in addition_combos:

                    nc_br = np.array(c_br,dtype=complex)

                    if len(cad) > 0:
                       
                        nc_br = replace_column(nc_br, cad, -1*col)

                    if check_is_diagram(nc_br, nmint, -1*col):
                        
                        for crm in removal_combos:

                            nc_ubr = np.array(c_ubr,dtype=complex)

                            if len(crm) > 0:

                                nc_ubr = replace_column(nc_ubr,crm, -1j*col)

                            if check_is_diagram(nc_ubr, t, -1j*col):
                                

                                if col == np.min(cols_br):
                                    
                                    n0, checking_diag,is_word = check_words(nc_br,nc_ubr,min_n0,return_n0=True)

                                    if is_word and (not str(checking_diag) in checking_diags):
                                        cands_br.append(nc_br)
                                        cands_ubr.append(nc_ubr)
                                        checking_diags.append(str(checking_diag))
                                        
                                else:
                                    cands_br.append(nc_br)
                                    cands_ubr.append(nc_ubr)
                                    
        candidates_br = np.array(cands_br)
        candidates_ubr = np.array(cands_ubr)
    

    cols_ubr = np.array(list(set(ubrA.real.astype(int).flatten())))
    cols_ubr = cols_ubr[cols_ubr>0]
    n0s = []

    if len(cols_ubr)==0:
        steps_up = np.zeros(len(candidates_br))
        n0s = np.ones(len(candidates_br))*min_n0
    else:
        checking_diags = []
        
    diag_labels = []
    
    for col in cols_ubr:
        
        num_boxes = np.sum(ubrA.real.astype(np.int64)==col)

        cands_br = []
        cands_ubr = []

        for cind in range(len(candidates_br)):
            c_br = candidates_br[cind]
            c_ubr = candidates_ubr[cind]
            tmax = min([np.sum(c_br!=0),num_boxes])

            for t in range(tmax+1):

                nmint = num_boxes - t
                where_to_add = get_entries_to_add(c_ubr,nmint)
                where_to_remove = get_entries_to_remove(c_br,t)

                if len(where_to_remove)>0:
                    
                    where_to_remove = where_to_remove[np.where(where_to_remove[:,1]<len(brB[0]))[0]]
                
                addition_combos = in_different_rows(np.array(list(it.combinations(where_to_add, nmint))))
                removal_combos = in_different_rows(np.array(list(it.combinations(where_to_remove, t))))
                if col == 1:
                    removal_combos = remove_combos_with_insufficient_gaps(removal_combos,c_br)
            

                for cad in addition_combos:

                    nc_ubr = np.array(c_ubr,dtype=complex)

                    if len(cad) > 0:

                        nc_ubr = add_column(nc_ubr, cad, 1*col)


                    if check_is_diagram(nc_ubr, nmint, col):
                        
                        for crm in removal_combos:

                            nc_br = np.array(c_br,dtype=complex)

                            n0_setp_up = 0
                            if len(crm) > 0:

                                nc_br = add_column(nc_br,crm, 1j*col)
                                

                            if check_is_diagram(nc_br, t, 1j*col):

                                
                                if col == max(cols_ubr):
                                    
                                    n0, checking_diag,is_word = check_words(nc_br,nc_ubr,min_n0,return_n0=True)
                                    
                                    if is_word and (not str(checking_diag) in checking_diags):
                                        n0s.append(max(n0,min_n0))
                                        cands_br.append(nc_br)
                                        cands_ubr.append(nc_ubr)
                                        checking_diags.append(str(checking_diag))
                                        diag_labels.append(checking_diag)

                                
                                elif check_words(nc_br,nc_ubr,min_n0):
                                    
                                    cands_br.append(nc_br)
                                    cands_ubr.append(nc_ubr)

        candidates_br = np.array(cands_br)
        candidates_ubr = np.array(cands_ubr)
        
        
    return unique_permutations_with_N0s(candidates_br,candidates_ubr,np.array(n0s))

