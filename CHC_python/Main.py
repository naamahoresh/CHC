import sys
import numpy as np
import random
from CHC_python import Pool
from CHC_python import Gene
from CHC_python import CHC
from CHC_python import read_config
from CHC_python import binary_function
from CHC_python import selecation_functions
from CHC_python import Status

import logging
import os
logging.basicConfig(level=logging.INFO, filemode='w', filename=os.getcwd()+'/logger.log', format='%(message)s')


'''
TODO:
1) understand if I need the dependency between the gene's length string and the number of params and number of bites to each param (gene's length string= # params * # bites represent param)
2) implement the fitness function
3) write the status better, like the old format
'''


StatusInterval = 100
DumpInterval = 10000
optimum_threshold = 9
MutateRate = 0.35
opt_operator="MAX"

def set_global_params():
def set_global_params():
    pass

def test_fun(numbers_array):
    numbers_array_tmp=np.copy(numbers_array)
    E=[]
    n=len(numbers_array_tmp)
    numbers_array_tmp[numbers_array_tmp == 0] = -1
    for k in range(1, n):
        E.append((numbers_array_tmp[0:n - k].dot(numbers_array_tmp[k:])) ** 2)

    fitness= (float(n ** 2) / float(2 * sum(E)))
    return fitness

if __name__ == "__main__":
    user_args=sys.argv
    pool_size=int(user_args[1])
    string_length = int(user_args[2])
    num_exp  =int(user_args[3])
    number_trials =int(user_args[4])
    selection_bias = float(user_args[5])

    opt_operator = "MAX"
    StatusInterval = 10
    DumpInterval = 10
    optimum_threshold = 9
    MutateRate = 0.35
    random.seed()
    final_results=[]

    param_file=os.getcwd()+"/F1/f1"
    result_file = os.getcwd()+"/F1/result"
    final_result_file = os.getcwd()+"/F1/final_result"

    # Coding, Cplxty, num_params, work_params = read_config.read_config(param_file)
    #
    # Status.save_params(work_params,num_params,result_file)

    #Main loop of the experiment
    for exp in range (num_exp):
        Status.print_progress_headers_by_exp (result_file, exp)

        if (exp==0):
            current_generation=0

        parent_pool = Pool.pool(pool_size, string_length)
        num_init_parent_pool = parent_pool.init_pool("",0,pool_size,test_fun)
        parent_pool.sort()

        org_threshold = threshold = string_length / 4
        eval_count = 0

        dict_best_result = {"worth":0, "num_trails":0 , "string":""}
        best_sol=0
        best_trail=0
        # Optimisation
        for trail in range (0,number_trials):
            children_produced = 0
            child_pool = Pool.pool(pool_size, string_length, 0)

            if (StatusInterval and trail%StatusInterval==0):
                Status.show_progress_exp(parent_pool,trail,result_file)

            #Combination loop
            for combination_index in range(int(pool_size/2)):

                mom_tmp, dad_tmp = parent_pool.get_parents(selection_bias, selecation_functions.linear)
                mom = Gene.gene(other_gene=mom_tmp)
                dad = Gene.gene(other_gene=dad_tmp)

                if (int(binary_function.hamming_distance (mom, dad, string_length)/2)>threshold):
                    CHC.hux(mom, dad,string_length)
                    mom.gene_worth = test_fun(mom.string)
                    dad.gene_worth = test_fun(dad.string)

                    child_pool.insert_gene(mom)
                    child_pool.insert_gene(dad)

                    children_produced += 2
                    eval_count += 2


            #selecation - choose the best from the parent and children
            num_child_insert = 0
            for insert_index in range (children_produced):
                if (selecation_functions.find_optimum(child_pool.genes_array[insert_index].gene_worth , parent_pool.genes_array[pool_size-1].gene_worth,opt_operator)):
                    num_child_insert+=1
                    parent_pool.replace_new_gene_by_worth_val(child_pool.genes_array[insert_index],opt_operator)

            if selecation_functions.find_optimum(parent_pool.genes_array[0].gene_worth, optimum_threshold,opt_operator):
                logging.info("The run found the optimum {} after {} iterations".format(parent_pool.genes_array[0].gene_worth,trail+1))
                break

            #restart the pop
            if (num_child_insert==0):
                threshold-=1
                if threshold<=0:
                    CHC.cataclysm(parent_pool, pool_size-1, string_length, 0, MutateRate,trail)
                    threshold = org_threshold

                    for index in range(pool_size):
                        # parent_pool.genes_array[index].worth = CHC.mstr_eval()  #m_eval
                        parent_pool.genes_array[index].gene_worth = test_fun(parent_pool.genes_array[index].string)  #m_eval
                        parent_pool.sort()
                        eval_count+=1

            if selecation_functions.find_optimum(parent_pool.genes_array[0].gene_worth , best_sol,opt_operator):
                best_sol=parent_pool.genes_array[0].gene_worth
                best_trail=trail

            # if (eval_count >= number_trials):
            #     Status.show_progress_exp(parent_pool,trail,result_file)
            #     break

            if (DumpInterval and (current_generation % DumpInterval==0)):
                Status.status_dump(parent_pool, param_file,user_args,trail)

        Status.show_progress_exp(parent_pool,trail,result_file)

        dict_best_result["worth"]=best_sol
        dict_best_result["num_trails"]=best_trail
        dict_best_result["string"]=str(parent_pool.genes_array[0].string)
        final_results.append(dict_best_result)


    Status.final_report(final_results,final_result_file) #stats


