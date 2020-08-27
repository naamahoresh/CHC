import math
import random


def linear(max, bias):
    if bias <= 1.0:
        index = int(random.uniform(0, 1) * max)
    else:
        index = int(
            max
            * (bias - math.sqrt(bias * bias - 4.0 * (bias - 1) * random.uniform(0, 1)))
            / 2.0
            / (bias - 1)
        )

    return index


def find_optimum(object_one, object_two, opt_operator):
    if opt_operator == "MAX":
        if object_one > object_two:
            return 1
    if opt_operator == "MIN":
        if object_one < object_two:
            return 1
    return 0
