cmake_minimum_required(VERSION 3.13)
project(CHC_CSUga C)

set(CMAKE_C_STANDARD 11)

include_directories(CHC/Evals)
include_directories(CHC/Include)
include_directories(Ga/Include)

add_executable(CHC_CSUga
        CHC/Evals/cclock.c
        CHC/Evals/ch.c
        CHC/Evals/cpi.c
        CHC/Evals/dejong_params.h
        CHC/Evals/f2c.h
        CHC/Evals/gasdev.c
        CHC/Evals/m_evals.c
        CHC/Evals/m_evals.h
        CHC/Evals/MOCMAParams.h
        CHC/Evals/ran2.c
        CHC/Evals/sleep.c
        CHC/Include/binarys.h
        CHC/Include/binarys_params.h
        CHC/Include/op_cataclysm.h
        CHC/Include/op_hux.h
        CHC/Include/op_hux_params.h
        CHC/Include/stats.h
        CHC/Include/stats_params.h
        CHC/Src/binarys.c
        CHC/Src/op_cataclysm.c
        CHC/Src/op_hux.c
        CHC/Src/stats.c
        CHC/main.c
        Ga/Include/ga_copy.h
        Ga/Include/ga_global.h
        Ga/Include/ga_global_extern.h
        Ga/Include/ga_params.h
        Ga/Include/ga_pool.h
        Ga/Include/ga_random.h
        Ga/Include/ga_selection.h
        Ga/Include/ga_signals.h
        Ga/Include/ga_status.h
        Ga/Include/ga_tournament.h
        Ga/Include/ga_xdr.h
        Ga/Include/gene.h
        Ga/Src/ga_copy.c
        Ga/Src/ga_params.c
        Ga/Src/ga_pool.c
        Ga/Src/ga_random.c
        Ga/Src/ga_selection.c
        Ga/Src/ga_signals.c
        Ga/Src/ga_status.c
        Ga/Src/ga_tournament.c
        Ga/Src/ga_xdr.c)
