import numpy as np
cimport numpy as cnp
from cython.cimports.libc.math import ceil
from collections import Counter

cdef inline bint is_monotonic_decreasing(cnp.ndarray[cnp.int64_t, ndim=1] lis):
    """Check if a list is monotonically decreasing."""
    cdef Py_ssize_t i
    cdef Py_ssize_t n = lis.shape[0]
    if n > 0:
        for i in range(n - 1):
            if lis[i] < lis[i + 1]:
                return False
    return True


def partition_tuplify(part):
    """Convert input into a tuple representation of a partition."""
    if not isinstance(part, tuple):
        if part == 0:
            return ()
        part = (part,)
    return part


def make_2d_array(perm, column_sort=False, dims=None):
    """
    Create a 2D array representation of a tableau from a partition.
    """
    cdef:
        cnp.ndarray[cnp.complex128_t, ndim=2] arr
        Py_ssize_t rows, cols, row, perm_row

    dim = (1, 1)
    if dims is not None:
        dim = dims

    arr = np.zeros(dim, dtype=np.complex128)
    if column_sort:
        arr *= (0.5 + 0.5j)

    if len(perm) > 0 or np.any(dim) > 1:
        rows = max(dim[0], len(perm))
        cols = max(dim[1], perm[0] if len(perm) > 0 else 0)

        arr = np.zeros((rows, cols), dtype=np.complex128)

        for row in range(len(perm)):
            perm_row = perm[row]
            arr[row, perm_row:] = 0  # Fill with 0s beyond the partition size
            if column_sort:
                arr[row, :perm_row] = np.arange(1, perm_row + 1, dtype=np.float64)
            else:
                arr[row, :perm_row] = np.ones(perm_row, dtype=np.complex128) * (0.5 + 0.5j)

    return np.array(arr)
    
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
    
def get_full_diag(cnp.ndarray[cnp.complex128_t, ndim=2] br, 
                  cnp.ndarray[cnp.complex128_t, ndim=2] ubr, 
                  int arg_min):
    """
    Concatenate the flipped portion of 'br' with the portion of 'ubr' up to 'arg_min'.
    """
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] flipped_br, concd

    # Flip 'br' twice: first along axis 1, then along axis 0
    flipped_br = np.flip(np.flip(br[:arg_min], axis=1), axis=0)

    # Concatenate flipped 'br' and 'ubr' along axis 1
    concd = np.concatenate((flipped_br, ubr[:arg_min]), axis=1)

    return concd
    
def col_add(tuple brA, tuple brB):
    """
    Add the first elements of two lists if they exist.
    """
    cdef int col = 0
    if len(brA) > 0:
        col += brA[0]
    if len(brB) > 0:
        col += brB[0]
    return col


def make_candidate_array(tuple pairB, tuple pairA):
    """
    Create candidate arrays from PairTableaux.
    """
    cdef:
        tuple brB, brA, ubrB, ubrA
        int rows_br, cols_br, rows_ubr, cols_ubr, cols
        cnp.ndarray[cnp.complex128_t, ndim=2] arr_ubr, arr_br

    # Extract barred and unbarred parts
    brB = pairB[0]
    brA = pairA[0]
    ubrB = pairB[1]
    ubrA = pairA[1]

    # Compute rows and columns
    rows_br = len(brB) + len(brA)
    cols_br = col_add(brA, brB)
    rows_ubr = len(ubrB) + len(ubrA)
    cols_ubr = col_add(ubrA, ubrB)

    # Determine the maximum number of columns
    cols = max(cols_br, cols_ubr)

    # Create the unbarred and barred arrays
    arr_ubr = make_2d_array(ubrB, dims=(rows_ubr + rows_br, cols), column_sort=False)
    arr_br = make_2d_array(brB, dims=(rows_br + rows_ubr, cols), column_sort=False) * -1

    return arr_br, arr_ubr


def make_pair_array(tuple pair, bint column_sort=False):
    """
    Create arrays for a pair with optional column sorting.
    """
    cdef:
        list part_br, part_ubr
        cnp.ndarray[cnp.complex128_t, ndim=2] arr_br, arr_ubr

    part_br = pair[0]
    part_ubr = pair[1]

    # Create the arrays
    arr_ubr = make_2d_array(part_ubr, column_sort=column_sort)
    arr_br = make_2d_array(part_br, column_sort=column_sort)

    return arr_br * -1j, arr_ubr

def is_subword(cnp.ndarray[cnp.float64_t, ndim=1] lis):
    """
    Check if a 1D complex array forms a subword.
    """
    cdef:
        list counter_lis, counter_keys
        int max_key
    
    counter_lis = list(Counter(lis).values())
    counter_keys = list(Counter(lis).keys())
    max_key = int(max(counter_keys)) if counter_keys else 0

    return is_monotonic_decreasing(np.asarray(counter_lis).astype(np.int64)) and \
        (sorted(counter_keys) == list(range(1, max_key + 1)))


def is_word(cnp.ndarray[cnp.float64_t, ndim=1] lis):
    """
    Check if a 1D complex array forms a word.
    """
    cdef:
        Py_ssize_t indx
        list checker

    if len(lis) > 1:
        checker = [is_subword(lis[0:indx + 1]) for indx in range(len(lis))]
        return np.all(checker)
    elif len(lis) == 1:
        return max(lis) <= 1
    else:
        return True


def check_wording(cnp.ndarray[cnp.complex128_t, ndim=2] concd):
    """
    Check the wording of a 2D complex array.
    """
    cdef:
        cnp.ndarray[cnp.complex128_t, ndim=2] br_word_array, ubr_word_array
        #np.ndarray[np.bool_, ndim=2] mask_br, mask_ubr
        cnp.ndarray[cnp.complex128_t, ndim=1] br_w, ubr_w
        cnp.ndarray[cnp.float64_t, ndim=1] br_word, ubr_word
        bint rtn

    # Flip and process br_word_array
    br_word_array = np.flip(concd, axis=0)
    mask_br = (br_word_array.real.astype(np.int64) < 0) | (br_word_array.imag.astype(np.int64) < 0)
    br_w = br_word_array[mask_br]

    br_word = np.zeros(len(br_w), dtype=np.float64)
    br_word[br_w.real < 0] -= br_w.real[br_w.real < 0]
    br_word[br_w.imag < 0] -= br_w.imag[br_w.imag < 0]

    # Flip and process ubr_word_array
    ubr_word_array = np.flip(concd.T, axis=1)
    mask_ubr = (ubr_word_array.imag.astype(np.int64) > 0) | (ubr_word_array.real.astype(np.int64) > 0)
    ubr_w = ubr_word_array[mask_ubr]

    ubr_word = np.zeros(len(ubr_w), dtype=np.float64)
    ubr_word[ubr_w.imag > 0] += ubr_w.imag[ubr_w.imag > 0]
    ubr_word[ubr_w.real > 0] += ubr_w.real[ubr_w.real > 0]

    # Combine results
    rtn = is_word(br_word) and is_word(ubr_word)
    return rtn

def check_shape(cnp.ndarray[cnp.complex128_t, ndim=2] diag):
    """
    Check if a tableau has a valid shape.
    """
    cdef:
        #cdef np.ndarray[np.bool_, ndim=2] mask
        cnp.ndarray[cnp.int64_t, ndim=2] dig
        cnp.ndarray[cnp.int64_t, ndim=1] perm
        Py_ssize_t rows, cols
        bint gaps
    
    rows = <Py_ssize_t>diag.shape[0]
    cols = <Py_ssize_t>diag.shape[1]
    mask = (diag.real != 0) & (diag.imag.astype(np.int64) <= 0)
    dig = <cnp.ndarray[cnp.int64_t, ndim=2]>np.empty((rows,cols), dtype=np.int64)
    dig[:, :] = 0  # Initialize to zero if required
    dig[mask] += 1

    perm = np.trim_zeros(np.sum(dig, axis=-1).astype(np.int64))
    if len(perm) > 0:
        gaps = True
        for col in range(dig.shape[1]):
            if not is_monotonic_decreasing(dig[:, col]):
                gaps = False
        #gaps = np.all(np.apply_along_axis(is_monotonic_decreasing, 0, dig))
        return is_monotonic_decreasing(perm) and gaps
    else:
        return True
        
def remove_combos_with_insufficient_gaps(removal_combos,diag):
    
    gaps_br = np.sum((diag.real == -1),axis =-1)
    
    recombs = np.array([removal_combos[comb] for comb in range(len(removal_combos))
        if np.all([(np.sum(gaps_br[0:removal_combos[comb,ind,0]+1]==0) >= len(removal_combos[comb,0:ind+1,0]))
                   for ind in range(len(removal_combos[comb]))]) 
                       ], dtype=np.int64)
    
    if len(recombs)==0:
        return np.array([[]], dtype=np.int64)
    
    return recombs


def in_different_rows(cands):
    """
    Filter candidates ensuring entries are in different rows.
    """

    try:
        cands_kept = np.array([cand for cand in cands if len(list(set(cand[:,0])))==len(cands[0,:,0])
                      ])
    
        return cands_kept
    except Exception as e:
        return np.array([[]])



