import logging
import os
import statistics
logging.basicConfig(level=logging.INFO, filemode='w', filename=os.getcwd()+'/logger.log', format='%(message)s')
from collections import Counter


def print_progress_headers_by_exp(result_file,exp_num):
    output_string = "\nExp num: {}\nGenerstion     Best        Worst       Median      Average     \n\n".format(exp_num)
    file = open(result_file,'a')
    file.write(output_string)
    file.close()

def show_progress_exp(pool, current_generation,result_file):
    avr_worth = 0
    for gene_index in range(pool.pool_size):
        avr_worth+=pool.genes_array[gene_index].gene_worth
    avr_worth=avr_worth/pool.pool_size

    output_string = "{:3d}            {:.4f}      {:.4f}      {:.4f}      {:.4f}     \n".format(
        current_generation,pool.genes_array[0].gene_worth,pool.genes_array[int(pool.pool_size-1)].gene_worth,pool.genes_array[int(pool.pool_size/2)].gene_worth,avr_worth)

    file = open(result_file,'a')
    file.write(output_string)
    file.close()

def status_dump(pool, file_name,user_args,trail):
    file_name_config = file_name+".config"
    file_name_pool = file_name+".pool"
    config_name = ["pool_size","string_length","num_exp","number_trials"]
    try:
        file_config = open(file_name_config, 'w')
    except:
        logging.error("Cannot open file to save the pool config: {}".format(file_config))
        exit(1)

    for config in range(1,len(user_args)):
            if config<len(config_name):
                file_config.write("{}: {}\n".format(config_name[config],user_args[config]))
            else:
                file_config.write("Unknown config: {}\n".format(user_args[config]))

    file_config.close()

    try:
        file_pool = open(file_name_pool, "w")
    except:
        logging.error("Cannot open file to save the pool: {}".format(file_pool))
        exit(1)

    poll_to_save = save_pool(pool,trail)
    file_pool.write(poll_to_save)
    file_pool.close()

def save_params(work_params,num_params,result_file):
    file = open(result_file,'w')
    file.write("Params:\n")
    for param in range(num_params):
        file.write(("Param Name: {} -- {}\n".format(work_params[param]['param_name'],work_params[param])))
    file.write("\n\n")
    file.close()


def save_pool(pool,trail):
    string ="Trail num: {}\n".format(trail)
    for gene_index in range(pool.pool_size):
        string_gene=""
        for string_index in pool.genes_array[gene_index].string:
            string_gene=string_gene+str(string_index)
        string = string+"Gene num: {}    Worth: {}    String: {}\n".format(str(gene_index),str((pool.genes_array[gene_index].gene_worth)),string_gene)
    return string

def final_report(dict_best_result, file_name):
    file = open(file_name,'w')
    worth_list=[]
    trail_list=[]
    string_list=[]
    for sol in dict_best_result:
        worth_list.append(sol["worth"])
        trail_list.append(sol["num_trails"])
        string_list.append(sol["string"])

    str_final = "Max solution: {} ({} exp found this result)\nMax String: {}\nAvarge: {}\nMedian: {}".format(max(worth_list), Counter(worth_list)[max(worth_list)],
                                                                                                            string_list[worth_list.index(max(worth_list))], statistics.mean(worth_list),
                                                                                                            statistics.median(worth_list))
    file.write(str_final)
    file.close()
    print(str_final)
