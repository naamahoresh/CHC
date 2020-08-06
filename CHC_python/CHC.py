import numpy as np
import math
import random
from CHC_python import binary_function
import logging
import os
logging.basicConfig(level=logging.INFO, filemode='w', filename=os.getcwd()+'/logger.log', format='%(message)s')

#TO-DO:
#To test: cataclysm, parm_cataclysm,ham_dist,gray, degray,transform,hux
#To implement: random_function and fracrand - implement in another calss. the Stats Functions



def hux(mom, dad, length):
    dif_array=[]
    count_diff=0
    for index in range (length):
        if mom.string[index]!=dad.string[index]:
            dif_array.append(index)
            count_diff+=1

    loop_size = round_int(float(count_diff/2))
    for index in range(loop_size):
        position= int(random.uniform(0, 1) * count_diff)
        while dif_array[position]<0:
            position = int(random.uniform(0, 1) * count_diff)

        indx = dif_array[position]
        temp = mom.string[indx]
        mom.string[indx] = dad.string[indx]
        dad.string[indx] = temp
        dif_array[position] = -1



def round_int(intager):
    intager_floor = math.floor(intager)
    if (intager - intager_floor >= .5):
        return (intager_floor + 1)
    else:
        return (intager_floor)


def cataclysm (pool, pool_size, string_length, seed_idx, pct_mut,trail):
    logging.info("Trail {}: We are inside the cataclysm function".format(trail))
    to_seed = seed_idx + 1

    for index_pop in range (pool_size):
        for index_string in range(string_length):
            is_mutate = int(random.uniform(0, 1)+pct_mut)

            if is_mutate:
                if pool.genes_array[seed_idx].string[index_string]==1:
                    pool.genes_array[to_seed].string[index_string] = 0
                else:
                    pool.genes_array[to_seed].string[index_string] = 1

            else:
                pool.genes_array[to_seed].string[index_string] = pool.genes_array[seed_idx].string[index_string]

        to_seed+=1

def parm_cataclysm(pool, size, nb_parms, sparse, perm_length, seed_idx, pct_mut):
    seed_prob =  sparse / perm_length
    to_seed = seed_idx + 1

    for m in range (size):
        index=0
        for k in range(nb_parms):
            mutate = random.uniform(0, 1)+pct_mut

            if mutate>0.5:
                for j in range (index, perm_length + index):
                    if random.uniform(0, 1)() <=  seed_prob:
                        pool[to_seed] = 1
                    else:
                        pool[to_seed] = 0

            else:
                for j in range (index, perm_length + index):
                    pool[to_seed] = pool[seed_idx]

        to_seed+=1





# def mstr_eval(Coding, Cplxty,num_params, work_params, string, str_length):#mstr_eval
#     if Coding=="gray":
#         new_str=binary_function.degray()
#         pass
#     else:
#         new_str=string
#
#     binary_function.sym_trnsfrm(new_str, num_params, work_params) #transform the binary to input parameters