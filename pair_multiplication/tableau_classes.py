import numpy as np
from .cutils_pair_multiplication import *
import itertools as it
from multiprocessing import Pool


def get_pair_base_labeled(arr):
    cand_r = arr.real#.astype(str)
    cand_i = arr.imag.astype(int)
    
    cand_r_pos = np.array(cand_r)
    cand_r_pos[cand_r.astype(int) < 0] = 0
    
    mask_none = (cand_r==0.0) + ((cand_i != 0)*(~(cand_r.astype(int)>0)))
    mask_empty = (cand_r==0.5) + (cand_r.astype(int)>0)
    mask_bullet = (cand_r==-0.5) + (cand_i < 0)

    string_arr = np.array(cand_r_pos).astype(int).astype(str)
    string_arr[mask_none] = r'\none'
    string_arr[mask_empty] = r'\empty'
    string_arr[mask_bullet] = r'\bullet'
    
    return string_arr


    

def get_barred_part_marked(arr):
    
    cand_i = arr.imag.astype(int)
    cand_r = arr.real.astype(int)
    imag_str_arr = cand_i.astype(str)
    
    negatives_re = cand_r < 0
    
    imag_str_arr[negatives_re] = np.char.add(r'\overline{', (-cand_r).astype(str)[negatives_re])
    imag_str_arr[negatives_re] = np.char.add(imag_str_arr[negatives_re], '}')

    imag_str_arr[imag_str_arr=='0'] = ''
    imag_str_arr[imag_str_arr!='0'] = np.char.add(imag_str_arr[imag_str_arr!='0'], '')
        
    negatives_im = (np.char.find(imag_str_arr, '-') == 0) * (~negatives_re)
    imag_str_arr[negatives_im] = np.char.add(r'\overline{', (-cand_i).astype(str)[negatives_im])
    imag_str_arr[negatives_im] = np.char.add(imag_str_arr[negatives_im], '}i')

    mask = negatives_im + negatives_re
    
    return imag_str_arr, mask



def get_unbarred_part_marked(arr):
    
    cand_i = arr.imag.astype(int)
    cand_r = arr.real.astype(int)
    
    imag_str_arr = cand_r.astype(str)
    
    positives_re = cand_r > 0
    
    imag_str_arr[positives_re] = cand_r.astype(str)[positives_re]
        
    positives_im = (cand_i > 0) * (~positives_re)
    imag_str_arr[positives_im] = cand_i.astype(str)[positives_im]
    imag_str_arr[positives_im] = np.char.add(imag_str_arr[positives_im], 'i')

    mask = positives_im + positives_re
    
    return imag_str_arr, mask



def process_iterate_cands_barred(args):
    
    col, num_boxes, min_n0, cand_pair, min_cols_br = args

    cands = CandidateList(min_n0=min_n0)

    tmax = min([np.sum(cand_pair.ubr.real!=0),num_boxes])
            
    for t in range(tmax+1):

        nmint = num_boxes - t
        addition_combos = cand_pair.br.get_entries_to_add(nmint)
        removal_combos = cand_pair.ubr.get_entries_to_remove(t)

        for cad in addition_combos:

            nc_br = cand_pair.br

            if len(cad) > 0:
                       
                nc_br = cand_pair.br.replace_column(cad, -1*col)

            if nc_br.check_is_diagram(nmint, -1*col):
                for crm in removal_combos:

                    nc_ubr = cand_pair.ubr

                    if len(crm) > 0:

                        nc_ubr = cand_pair.ubr.replace_column(crm, -1j*col)

                    if nc_ubr.check_is_diagram(t, -1j*col):
                        if col == min_cols_br:
                            
                            pair = PairTableau(tableau=[nc_br,nc_ubr], complete = True, min_n0=min_n0)

                            if pair.admissible:
                                
                                cands.append(pair)
                                        
                        else:
                            pair = PairTableau(tableau=[nc_br,nc_ubr], complete = False, min_n0=min_n0)
                            
                            if pair.admissible:
                                
                                cands.append(pair)
                            
                            
    return cands

def process_iterate_cands_unbarred(args):
    
    col, num_boxes, min_n0, cand_pair, max_cols_ubr, len_brb = args
    
    cands = CandidateList(min_n0=min_n0)
    
    tmax = min([np.sum(cand_pair.br.array!=0),num_boxes])

    for t in range(tmax+1):

        nmint = num_boxes - t
        addition_combos = cand_pair.ubr.get_entries_to_add(nmint)
        removal_combos = cand_pair.br.get_entries_to_remove(t, len_brb = len_brb)

        if col == 1:
            removal_combos = remove_combos_with_insufficient_gaps(removal_combos,cand_pair.br.array)

        for cad in addition_combos:

            nc_ubr = cand_pair.ubr

            if len(cad) > 0:
                nc_ubr = cand_pair.ubr.add_column(cad, 1*col)

            if nc_ubr.check_is_diagram(nmint, col):
                        
                for crm in removal_combos:

                    nc_br = cand_pair.br

                    if len(crm) > 0:
                        nc_br = cand_pair.br.add_column(crm, 1j*col)
                                
                    if nc_br.check_is_diagram(t, 1j*col):
      
                        if col == max_cols_ubr:
                                    
                            pair = PairTableau(tableau=[nc_br,nc_ubr], complete = True, min_n0=min_n0)

                            if pair.admissible:
                                
                                cands.append(pair)

                                
                        else:
                            pair = PairTableau(tableau=[nc_br,nc_ubr], complete = True,min_n0=min_n0)
                            
                            if pair.admissible:
                                
                                cands.append(pair)
                            
    return cands


def parallel_process_iterate_cands_barred(candidates, num_boxes, min_n0, cols_br, col):
    
    min_cols_br = np.min(cols_br)
    args = [[col, num_boxes, min_n0, candidates[ind], min_cols_br]
            for ind in range(candidates.length)]
    
    cands = CandidateList(min_n0=min_n0)
    
    with Pool(None) as pool:

        for result in pool.imap(process_iterate_cands_barred, args):
            
            cands += result
        
    return cands


def parallel_process_iterate_cands_unbarred(candidates, num_boxes, 
                                            min_n0, cols_ubr, col, len_brb):
    
    max_cols_ubr = np.max(cols_ubr)
    args = [[col, num_boxes, min_n0, candidates[ind], max_cols_ubr, len_brb]
            for ind in range(candidates.length)]
    
    cands = CandidateList(min_n0=min_n0)
    
    with Pool(None) as pool:

        for result in pool.imap(process_iterate_cands_unbarred, args):

            cands += result

    return cands




def pair_multipy(prB,prA):
    
    pairA = PairTableau(partitions = prA,column_sort=True,working_pair=False)
    pairB = PairTableau(partitions = prB, column_sort=True,working_pair=False)
    len_brb = len(pairB.br.array[0])
    
    min_n0 = max(len(prA[1])+len(prA[0]),len(prB[0])+len(prB[1]))
    
    candidates = CandidateList(partition_pair = [prB,prA], min_n0 = min_n0)

    cols_br = np.array(list(set(-pairA.br.imag.astype(int).flatten())))
    
    cols_br = np.flip(cols_br[cols_br>0])
    
    for col in cols_br:

        num_boxes = np.sum(-pairA.br.imag.astype(int)==col)

        candidates = parallel_process_iterate_cands_barred(candidates, 
                                                           num_boxes, min_n0, 
                                                           cols_br, col)
    
    cols_ubr = np.array(list(set(pairA.ubr.real.astype(int).flatten())))
    cols_ubr = cols_ubr[cols_ubr>0]
    
    for col in cols_ubr:
        
        num_boxes = np.sum(pairA.ubr.real.astype(np.int64)==col)

        candidates = parallel_process_iterate_cands_unbarred(candidates, 
                                                             num_boxes, min_n0, 
                                                             cols_ubr, col,len_brb)
    
    return candidates
    
    
    
def pair_multipy_verbosely(prB,prA):
    
    pairA = PairTableau(partitions = prA,column_sort=True,working_pair=False)
    pairB = PairTableau(partitions = prB, column_sort=True,working_pair=False)
    len_brb = len(pairB.br.array[0])
    
    min_n0 = max(len(prA[1])+len(prA[0]),len(prB[0])+len(prB[1]))
    
    candidates = CandidateList(partition_pair = [prB,prA], min_n0 = min_n0)

    cols_br = np.array(list(set(-pairA.br.imag.astype(int).flatten())))
    
    cols_br = np.flip(cols_br[cols_br>0])
    
    candidates_in_steps = []
    
    for col in cols_br:

        num_boxes = np.sum(-pairA.br.imag.astype(int)==col)

        candidates = parallel_process_iterate_cands_barred(candidates, 
                                                           num_boxes, min_n0, 
                                                           cols_br, col)
                                                           
        candidates_in_steps.append(candidates)
    
    cols_ubr = np.array(list(set(pairA.ubr.real.astype(int).flatten())))
    cols_ubr = cols_ubr[cols_ubr>0]
    
    for col in cols_ubr:
        
        num_boxes = np.sum(pairA.ubr.real.astype(np.int64)==col)

        candidates = parallel_process_iterate_cands_unbarred(candidates, 
                                                             num_boxes, min_n0, 
                                                             cols_ubr, col,len_brb)
        candidates_in_steps.append(candidates)
    return candidates_in_steps


class Tableau():
    
    def __init__(self, partition = None, tableau = None, barred = False, set_barred = False, column_sort=False):
        
        self.partition = None
        if not partition is None:
            self.partition = partition_tuplify(partition)
            
        self.array = None
        
        if tableau is None:
            self.array = make_2d_array(self.partition, column_sort = column_sort)
        else:
            self.array = np.array(tableau)
            
        if barred:
            
            if set_barred:
                self.array *= -1j
        
        self.barred = barred
            
        self.admissible_shape = check_shape(self.array) * check_shape(self.array.T)
        self.real = self.array.real
        self.imag = self.array.imag
        self._hash = hash((self.array.tobytes(), self.admissible_shape))
        
    def get_partition(self):

        mask_br = (self.real!=0) * (self.imag.astype(np.int64)<=0)
        shape_arr = np.zeros(self.array.shape)
        shape_arr[mask_br] += 1

        perm_arr = np.sum(shape_arr,axis = -1)
        
        perm = tuple(np.trim_zeros(perm_arr).astype(np.int64))
        
        self.permutation = perm
        
        return perm
            
    def check_is_diagram(self, num_entries, val):
        
        if self.admissible_shape:
    
            mask = (self.real!=0) * (self.imag.astype(np.int64)<=0)

            dig = np.zeros(self.array.shape, dtype=np.int64)
            dig[mask] += 1

            perm = np.trim_zeros(np.sum(dig,axis = -1).astype(np.int64))

            checker_num_added = False
            if val.real != 0:
                checker_num_added = np.sum(self.real == val.real) == num_entries
            elif val.imag != 0:
                checker_num_added = np.sum(self.imag == val.imag) == num_entries

            return checker_num_added
        else:
            return False
        
    def replace_column(self, adding_col, label):
    
        diag = np.array(self.array)
        rows,cols = adding_col.T
        diag[rows,cols] = label
        
        return Tableau(tableau = diag, barred = self.barred)


    def add_column(self, adding_col, label):

        diag = np.array(self.array)
        rows,cols = adding_col.T

        diag[rows,cols] = diag[rows,cols].real.astype(np.int64)+1j*diag[rows,cols].imag.astype(np.int64)

        diag[rows,cols] += label

        return Tableau(tableau = diag, barred = self.barred)
    
    
    def __hash__(self):
    
        return self._hash
    
    
    def get_entries_to_add(self,num):
    
        arr = self.array.copy()#np.array(array, dtype=np.complex128)
        coast_arr = arr.real==0
        coast = np.argwhere(coast_arr==True)

        if len(coast) > 0:

            ents = np.array([cl
                             for cl in coast 
                         if (cl[0] < np.min((coast[:,0])[coast[:,1]==cl[1]])+num)])
            if len(ents)>0:
                return in_different_rows(np.array(list(it.combinations(ents, num))))
            else:
                return np.array([[]])
        else:
            return np.array([[]])


    def get_entries_to_remove(self,num, len_brb = None):

        arr = self.array.copy()#np.array(array)
        coast_arr = (arr.imag.astype(np.int64)==0) * (arr.real!=0)
        coast = np.argwhere(coast_arr==True)
        
        if len(coast)>0:
            ents = np.array([cl
                             for cl in coast 
                         if (cl[0] >= np.max((coast[:,0])[coast[:,1]==cl[1]])-num+1)])
            if len(ents)>0:
                if not len_brb is None:
                    ents = ents[np.where(ents[:,1]<len_brb)[0]]

                return in_different_rows(np.array(list(it.combinations(ents,num))))
            else:
                return np.array([[]])
        else:
            return np.array([[]])


    
class PairTableau:
    
    def __init__(self, partitions = None, tableau = None, complete = True, min_n0 = 0, column_sort=False, 
                 working_pair = True):

        self.br = None
        self.ubr = None
        self.n0 = min_n0
        
        if not tableau is None:
            self.br = tableau[0]
            self.ubr = tableau[1]
        elif partitions is not None:
            self.br = Tableau(partition = partitions[0], barred = True, set_barred=True,
                              column_sort = column_sort)
            self.ubr = Tableau(partition = partitions[1], 
                               column_sort = column_sort)
            
        self.composite = None
        if working_pair:
            self.composite = self.get_overlap()
        else:
            complete = False
        
        self.admissible_word = True
        
        if complete:
            self.admissible_word = self.check_ones()

        self.complete = complete
        
        self.admissible_shape = self.br.admissible_shape * self.ubr.admissible_shape
        
        self.admissible = self.admissible_shape * self.admissible_word
        self._hash = hash((self.br,self.ubr, self.n0))
        
    def string_candidate(self):
    
        if not self.composite is None:
            n0 = self.n0
            
            comp = self.composite
            string_arr = get_pair_base_labeled(comp)
        
            barred, mask_br = get_barred_part_marked(comp)
            unbarred, mask_ubr = get_unbarred_part_marked(comp)
            
            num_rows = len(string_arr)
            label_str = ['']*num_rows
            
            append_col = True

            if not np.any((comp.real.astype(int) > 0) + (comp.imag.astype(int) > 0), axis = (0,1)):

                string_arr[mask_br] = barred[mask_br]
                append_col = False
            else:

                label_str = np.array([r'\none']*num_rows,dtype='<U21')
                plus = np.any(barred == r'\overline{1}i', axis = 1)
                label_str[plus] = np.char.add(label_str[plus],r'[\boxplus]')
                times = np.any(barred == r'\overline{1}', axis = 1)
                label_str[times] = np.char.add(label_str[times],r'[\boxtimes]')
                
                label_str[((~times)*(~plus))] = np.char.add(label_str[((~times)*(~plus))],r'[\circ]')
                
                string_arr[mask_ubr] = unbarred[mask_ubr]
                
            selector = np.any(string_arr!=r'\none',axis = 0)
            string_arr = ((string_arr.T)[selector]).T
            
            if append_col:
                string_arr = np.concatenate((label_str.reshape((num_rows,1)),string_arr),axis = 1)
                
            string_arr = np.array([[item.replace('0',r'\bullet') if r'\none[\boxtimes]' in row
                                    else 
                                    item.replace('0',r'\none')
                                    for item in row] for row in string_arr])
            
            full_strin = '\\\\\n'.join([r' & '.join(list(row))for row in string_arr])
                
            return r'\begin{ytableau}'+'\n'+full_strin+'\n'+r'\end{ytableau}'+'\n'
        else:
            return ''
            
    def __lt__(self,other):
    
        try:
            if self.n0 < other.n0:
                return True
            elif self.n0 == other.no:
                return len(self.string_candidate()) < len(other.string_candidate())
            else:
                return False
        except Exception as e:
            return False
            
            
    def __hash__(self):
        
        return self._hash
    

    def get_overlap(self):

        first_cols=np.array([self.br.array[:,0], self.ubr.array[:,0]])
        arg_min = np.sum(np.argmin(abs(first_cols),axis = 1),axis = 0)

        concd = np.concatenate((np.flip(np.flip(self.br.array[0:arg_min],axis = 1),
                                        axis =0),self.ubr.array[0:arg_min]), axis = 1)

        c_ubr = np.flip(np.trim_zeros(first_cols[1],'b'))

        c_br = np.flip(np.trim_zeros(first_cols[0],'b'))

        kmax = min(len(c_ubr),len(c_br))

        overlap = max([0]+[k+1 for k in range(kmax)
                if columns_align(np.flip(c_br[0:k+1]),c_ubr[0:k+1]) and
                check_wording(get_full_diag(self.br.array,self.ubr.array,arg_min-(k+1)))
                          ])

        arg_min -= overlap
        
        n0_base = len(self.br.get_partition())+len(self.ubr.get_partition())
        
        concd = get_full_diag(self.br.array,self.ubr.array,max(max(arg_min,n0_base),self.n0))

        self.n0 = len(concd.T)

        if not check_wording(concd):
            concd = None

        return concd
    
    
    def check_ones(self):

        rtn = False
        if not (self.composite is None):
            rtn = True

        if rtn:
            br_word_array = np.flip(self.composite,axis = 0)
            ubr_word_array = np.flip(self.composite.T,axis = 1)

            n0n = len(ubr_word_array[0])

            min1s = np.sum(ubr_word_array.real == -1,axis = 0)+\
                    np.sum(ubr_word_array.imag == -1,axis = 0)

            ones = np.sum(ubr_word_array.real == 1,axis = 0)+\
                   np.sum(ubr_word_array.imag == 1,axis = 0)

            step_up = [0]+[np.sum(ones[0:k+1])-np.sum(min1s[0:k+1]==0) for k in range(len(ones)) 
                              ]

            n0n += max(step_up)
            
            self.n0 = n0n

            return rtn
        else:
            return False

        
class CandidateList():
    
    def __init__(self,partition_pair = None, candidate_lis = None, min_n0 = 0):
        
        self.pairs = []
        self.n0s = []
        self.length = 0
        
        if not partition_pair is None:
            
            diagC_br, diagC_ubr = make_candidate_array(*partition_pair)
        
            pair = PairTableau(tableau = [Tableau(tableau = diagC_br, barred = True, set_barred = True),
                                          Tableau(tableau = diagC_ubr, barred = False)], min_n0=min_n0)
            self.pairs = [pair]
            self.n0s = [pair.n0]
            self.length = 1
            
        elif not candidate_lis is None:
            
            self.pairs = list(set(candidate_lis))
            self.length = len(self.pairs)
            self.n0s = [pair.n0 for pair in self.pairs]
            
            
    def __getitem__(self, index):
        
        return self.pairs[index]
    
    def append(self, candidate):
        
        if candidate not in self.pairs:
            
            self.pairs.append(candidate)
            self.n0s.append(candidate.n0)
            self.length = len(self.pairs)
        
    def __add__(self, cands_lis):
        
        new_cands = [cand for cand in cands_lis if (cand not in self.pairs)]
        self.pairs.extend(new_cands)
        
        self.n0s.extend([cand.n0 for cand in new_cands])
        self.length = len(self.pairs)
        
        return self
    
    def __radd__(self, other):
        
        return self+other
    
    
    def get_partitions(self):
        
        return [(self.pairs[ind].br.get_partition(), self.pairs[ind].ubr.get_partition()) 
                for ind in range(self.length)]
    
    def get_n0s(self):
        
        return self.n0s
        
    def string_list(self):
    
        joined = r'\,\oplus\,'.join([self.pairs[ind].string_candidate() 
                            for ind in range(self.length)]).replace(r'\,\oplus\,\,\oplus',r'\,\oplus')
        return joined.replace(r'\,\oplus\,', r'\,\oplus\, '+'\n')
    
