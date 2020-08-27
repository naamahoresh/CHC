import logging
import os

logging.basicConfig(
    level=logging.INFO,
    filemode="w",
    filename=os.getcwd() + "/logger.log",
    format="%(message)s",
)


def read_config(file_name):
    file_name = file_name + ".params"
    try:
        file = open(file_name, "r")
    except:
        logging.error("Cannot open parameter configuration file: {}".format(file_name))
        exit(1)

    function = file.readline()
    if not (function[0] == "f"):
        logging.error(
            "Need to implement function with noise, composition of functions and dual function evaluation"
        )
        exit(1)

    Coding = file.readline()[:-1]  # coding method
    if not (Coding == "bcd" or Coding == "gray"):
        logging.error("Unfamiliar coding type")
        exit(1)
    logging.info("Coding selected: {}".format(Coding))

    Cplxty = file.readline()[:-1]  # complexities
    if Cplxty != "simple":
        logging.error("Use only simple / standard complexities level")
        exit(1)

    Nb_params = int(file.readline()[:-1])
    logging.info("NUmber of params: {}".format(str(Nb_params)))

    postn = 0
    work_params = []
    for index in range(Nb_params):
        param = file.readline()
        params_list = param.split()
        if params_list[2] == "u":
            sign = 0
        else:
            sign = 1
        dic = {
            "param_name": params_list[0],
            "length": int(params_list[1]),
            "sign": sign,
            "scale": float(params_list[3]),
            "wrap_low": float(params_list[4]),
            "wrap_high": float(params_list[5]),
            "position": postn,
        }
        postn += dic["length"]
        dic["scale"] = 1 / dic["scale"]
        work_params.append(dic)
    return (Coding, Cplxty, Nb_params, work_params)


# read_config("/home/naamah/Documents/CHC/CHC_python/F1/f1.params")
