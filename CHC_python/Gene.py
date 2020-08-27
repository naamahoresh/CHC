import numpy as np
import random


class gene:
    def __init__(
        self,
        range_for_init=2,
        gene_length=0,
        translation_funcation=None,
        is_bias=0,
        bias=0.5,
        other_gene=None,
    ):

        if not other_gene == None:
            self.string = np.copy(other_gene.string)
            self.gene_worth = other_gene.gene_worth

        else:
            if not is_bias:
                self.string = np.random.randint(0, range_for_init, gene_length)
            else:
                self.string = np.zeros(gene_length)
                for index in range(gene_length):
                    if random.uniform(0, 1) < bias:
                        self.string[index] = 1

            self.gene_worth = translation_funcation(self.string)
