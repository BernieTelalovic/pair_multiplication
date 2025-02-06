# utils.py
import warnings
from collections import Counter
import numpy as np
import itertools as it
from .cutils_pair_multiplication import *
from multiprocessing import Pool

def Nc_positive(Nc:int):
        
    Nc = max([Nc,0])
        
    return Nc

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


def is_column_increasing(lists):
    
    return np.array([np.all([is_serial_increasing(li) for li in 
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
#########################################################################################################


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
        
        return np.array(arr)
    
    

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



def check_is_diagram(diag, num_entries, val):
    
    mask = (diag.real!=0) * (diag.imag.astype(np.int64)<=0)

    dig = np.zeros(diag.shape, dtype=np.int64)
    dig[mask] += 1
    
    perm = np.trim_zeros(np.sum(dig,axis = -1).astype(np.int64))
    
    checker_num_added = False
    if val.real != 0:
        checker_num_added = np.sum(diag.real == val.real) == num_entries
    elif val.imag != 0:
        checker_num_added = np.sum(diag.imag == val.imag) == num_entries
    
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
    
################################################ LaTeXing ##############################################################
########################################################################################################################


def fix_name(pair,name):
    
    if not name:
        name = 'you forgot to name me'
    
    return toCamel(name)
    
def toCamel(s):
    s = s.replace('_', ' ')
    ret = ''.join(x for x in s.title() if not x.isspace())     
    return ret[0].lower() + ret[1:]
    
def enumerate_YD(part):

    pair_part = ((),part)
    
    return enumerate_pair(pair_part)
    
def enumerate_YD_barred(part):

    pair_part = (part,())
    
    return enumerate_pair(pair_part)
    
def compose_ytab_from_array(arr):
    
    rows = []
    for row in arr:
        if not np.all(row==r'\none'):
            rw = '&'.join(row)
            rows.append(rw+r'\\')
    
    strin = '\n'.join(rows)
    
    if len(strin)>0:
        return r'\begin{ytableau}'+'\n'+strin+'\n'+r'\end{ytableau}'+'\n'
    else:
        return ''
    
def enumerate_pair(pair):

    nc_br = int(np.max(pair[0],initial=0))
    nc_ubr = int(np.max(pair[1],initial=0))

    nr_br = len(pair[0])
    nr_ubr = len(pair[1])

    nrows = nr_br+nr_ubr
    ncols = nc_br+nc_ubr

    barred_left = r'\overline{'
    barred_right = r'}'

    rev_ubr = np.array(pair[0])[::-1]

    rows = []
    for r in range(nrows):
        
        row = []
        if r < nr_ubr:
            
            nones = [r'\none']*nc_br
            labs = list(np.arange(1,pair[1][r]+1).astype('str'))
            more_nones = [r'\none']*(nc_ubr-pair[1][r])
            
            row = nones+labs+more_nones

            rows.append(row)
        else:
            labs = list(np.arange(-rev_ubr[r-nr_ubr],0))
            
            nones = [r'\none']*(nc_br-rev_ubr[r-nr_ubr])
            more_nones = [r'\none']*(nc_ubr)
            
            lab_strin = [barred_left+str(abs(lb))+barred_right for lb in labs]
            
            row = nones+lab_strin+more_nones

            rows.append(row)
    
    return np.array(rows)
    
    
def latex_pair(pair, name=None):
    
    name = fix_name(pair,name)
    
    ubr = pair[1]
    br = pair[0]
    
    num_ubr_rows = len(ubr)
    
    num_br_rows = 0
    if len(br) > 0:
        num_br_cols = max(br)
        
    br_str = '0,'*num_ubr_rows+','.join([f"{num_br_cols-row}+{row}" for row in list(reversed(br))])
    
    ubr_str = ','.join([f"{num_br_cols}+{row}" for row in ubr])
    
    print(r"\newcommand{"+"\\"+name+"}{\ydiagram"+r"[*(white)\bullet]{"+\
            br_str+"}*[*(white)]{"+ubr_str+"}}")
