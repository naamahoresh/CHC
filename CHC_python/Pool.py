import numpy as np
from CHC_python import Gene
from CHC_python import selecation_functions
import logging
import os
logging.basicConfig(level=logging.INFO, filemode='w', filename=os.getcwd()+'/logger.log', format='%(message)s')

# TODO:
# not implement yet: init_pool_from_file (if I will found an example file)
#
# Need to decide if:
# 1) the Gene will be an int or a string (right now, int)
# 2) gene be a separate class of a string\np array (right now, separate class)
# 3) init multiple types of population.(right now, only binary pop)


class pool:

    def __init__(self, pool_size_,string_length_,numer_genes_in_pool_=-1):
        self.genes_array = []
        self.pool_size = pool_size_
        self.string_length=string_length_
        if numer_genes_in_pool_>=0:
            self.numer_genes_in_pool = numer_genes_in_pool_
        else:
            self.numer_genes_in_pool = pool_size_



    def init_pool(self,init_data_file, start_index, stop_index, translation_function,bias=0):
        if (len(init_data_file)==0):
            if bias==0:
                num_init = self.init_random_pool(start_index, stop_index, translation_function)
            else:
                num_init = self.init_bias_pool(bias,start_index, stop_index, translation_function)

        else:
            try:
                file = open(init_data_file, "r")
            except:
                logging.ERROR("Can not open pool initialization file: {}".format(init_data_file))
                exit()

            num_init = self.init_pool_from_file(file, start_index, stop_index, translation_function)
            if num_init<0:
                logging.ERROR("Init pool: bad initialization indices")
                exit()

            if (num_init != stop_index-start_index): #If file did not contain enough initialization values, finish off with random initialization
                logging.WARNING("{} genes read from file '{}'.\n {} genes will be randomly generated\n".format(num_init, init_data_file, stop_index-num_init))
                rand_num_init = self.init_random_pool(start_index, stop_index, translation_function)
                num_init = num_init + rand_num_init

            file.close()

        return num_init

    """
    FUNCTION: init_random_pool for binary numbers
    DESCRIBE: Randomly initializes the array data structures of a genetic pool;
    INPUT PARAMETERS: start & stop position in pool;
                      pointer to evaluation function;
    RETURN VALUE: number of genes initialized
    """
    def init_random_pool(self, start_index, stop_index, translation_function):

        # Init Type:
        # binary = 2
        # scale = N

        init_type = 2

        for index in range (start_index, stop_index):
            tmp_gene=Gene.gene(init_type, self.string_length,translation_function) #if I want to init non binary population, I should change the 2
            self.genes_array.append(tmp_gene)

        return stop_index-start_index

    """
    FUNCTION: init_random_pool for binary numbers
    DESCRIBE: Initializes the array data structures of a binary genetic pool according to some biased percentage of 1's;
    INPUT PARAMETERS: bias fraction;
                      start & stop position in pool;
                      pointer to evaluation function;
    RETURN VALUE: number of genes initialized
    """

    def init_bias_pool(self, bias, start_index, stop_index, translation_function):
        # Init Type:
        # binary = 2
        # scale = N

        init_type = 2

        for index in range (start_index, stop_index):
            tmp_gene = Gene.gene(init_type, self.string_length,translation_function,1,bias)
            self.genes_array.append(tmp_gene)
        return stop_index-start_index



    def sort(self,opt="MAX"):
        worth_array = np.zeros(self.pool_size)
        for index in range(self.pool_size):
            worth_array[index]=self.genes_array[index].gene_worth

        sort_indexs = worth_array.argsort()
        string_tmp = []
        for index_worth in (range(len(self.genes_array))):
            string_tmp.append(self.genes_array[sort_indexs[index_worth]])

        if opt=="MAX":
            string_tmp=string_tmp[::-1]

        self.genes_array = string_tmp


    def init_pool_from_file(self,file, strt_index, stp_index, translation_function):

        print("Initialized from file isn't implement yet. Genes will generated randomly")
        return 0


    def get_parents(self, bias, bias_fun): #maybe to put it in the selection module??
        mom = bias_fun(self.pool_size, bias)
        dad = bias_fun(self.pool_size, bias)

        if (self.pool_size > 1):
            while (mom == dad):
                dad = bias_fun(self.pool_size,bias)

        return self.genes_array[mom],self.genes_array[dad]

    def insert_gene(self, new_gene):
        if len(new_gene.string)!=self.string_length:
            logging.error("The new gene doesn't have the same string length as the genes in the pool")
            exit(1)
        if len(self.genes_array)>= self.pool_size:
            logging.error("The pool is full (Poll size: {}), can't insert new gene".format(self.pool_size))
            exit(1)
        self.genes_array.append(new_gene)
        self.numer_genes_in_pool=+1

    def replace_new_gene_by_worth_val(self, new_gene,opt_operator):
        if not (selecation_functions.find_optimum(new_gene.gene_worth ,self.genes_array[self.pool_size-1].gene_worth,opt_operator)):
            return
        self.genes_array[self.pool_size-1]=new_gene
        self.sort()



def next_prime_number( n):
    next_perime=n+1
    while not is_prime(next_perime):
        next_perime+=1
    return next_perime

def is_prime( number):
    for index in range(2,int(number/2)):
        if number % index == 0:
            return 0
    else:
        return 1

def test_fun(numbers_array):
    return sum(numbers_array)



def pool_tester():
    #test the random_init_pool
    parent_pool = pool(5,10)
    parent_pool.init_pool("",0,5,test_fun)
    #
    for gene in parent_pool.genes_array:
        print(gene.gene_worth)
    print("\n")

    #test the sort
    parent_pool.sort()
    for gene in parent_pool.genes_array:
        print(gene.gene_worth)
    print("\n")

    #test the bias_init_pool
    parent_pool_bias = pool(5,10)
    parent_pool_bias.init_pool("",0,5,test_fun,0.8)

    for gene in parent_pool_bias.genes_array:
        print(gene.gene_worth)
    print("\n")

    #test prime number
    print(next_prime_number(11))
    print("\n")

    print(parent_pool.get_parents(bias=1.0, bias_fun=selecation_functions.linear))


# pool_tester()