import numpy as np
import math

def hamming_distance(Gene1, Gene2,array_length):
    distance = 0
    for index in range(array_length):
        if (Gene1.string[index] != Gene2.string[index]):
            distance +=1

    return distance


def hamming_diff(Gene1, Gene2,array_length):

    distance = 0
    for index in range(array_length):
        if (Gene1.string[index] != Gene2.string[index]):
            distance += 1

        if (distance>1): #HAM_LIMIT
            return 0

    return 1

def ham_dist(array1, array2,array_length,Work_params):
    distance = 0
    for index in range(array_length):
        delta1 = char_to_double(array1[Work_params[index].posn], Work_params[index].lngth,1,0) #NO_CHANGE, UNSIGNED
        delta2 = char_to_double(array2[Work_params[index].posn], Work_params[index].lngth,1,0) #NO_CHANGE, UNSIGNED

        distance+=abs(delta1 - delta2)
        if (distance>1): #HAM_LIMIT
            return 0

    return 1


#ctod
def char_to_double(string, length, step_value, sign):
    sum = 0
    sign_val=0

    if (((length < 2) and sign) or (length < 1)):
        return sum

    if sign:
        if string[0]=="1":
            sign_val=1

    for index in range (1,length): ##Naama - start from 1 or 0?
        sum+=sum
        if string[index]=="1":
            sum += 1

    if sign_val:
        sum = (-1) * sum

    if sign_val:
        sum = sum - 1

    sum = sum * step_value

    return sum


def gray(array_in, array_length,Work_params):
    array_out = np.zeros(len(array_in))
    for index in range (array_length):
        lngth = Work_params[index].lngth
        last = "0"
        for index2 in range (lngth):
            if array_in[index2]!=last:
                array_out[index2]=="1"
            else:
                array_out[index2]=="0"
            last= array_in[index2]
    return array_out


def degray(array_in, array_length,Work_params):
    array_out = np.zeros(len(array_in))
    for index in range (array_length):
        lngth = Work_params[index].lngth
        last = "0"
        for index2 in range (lngth):
            if array_in[index2]=="1":
                if last=="0":
                    last=="1"
                else:
                    last=="0"
            array_out[index2]=last
            last= array_in[index2]##check in the code. it write: last= array_out[index2]

    return array_out

def sym_trnsfrm(array, num_params,Work_params):
    values=np.zeros(num_params)
    buff_index=0
    for index in range (num_params):
        values[index] = char_to_double(array,Work_params[index].position, Work_params[index].length,
                                             Work_params[index].scale, Work_params[index].sign)

        buff_index+=Work_params[index].length
        max=math.pow(2, Work_params[index].length)
        values[index] = ((values[index] * Work_params[index].wrap_high) / (max / 2.0)) - Work_params[index].wrap_high

    return values
