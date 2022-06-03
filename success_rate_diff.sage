from scipy import stats
import numpy as np


log2 = lambda x: log(1.0*x)/log(2.0)
Phi = lambda x: (1 +erf(x/sqrt(2 )))/2 
Phi_inv = lambda x: sqrt(2 )*erfinv(2 *x-1 )

#  para for computing time is (keysize, K_g, l_0, l_r) 

distinguishers = {
    #'Asiacrypt21: Simeck64/128  r30 A': (64, 2^64, 2^2.59, 2^-1, 6, 2^60.41, [51, 128, 123, 9, 122] ,[]),

    #'Asiacrypt21: Simeck64/128  r30 B': (64, 2^64, 2^2.59, 2^-1, 5, 2^60.41, [51, 128, 123, 9, 122] ,[]),



    'Simon48/72  r17 A': (48, 2^47.5, 2^1.01, 2^-1.5, 1, 2^-45.49, [72, 69, 7+30, 41, 18]),
    #'Simon48/72  r17 B': (48, 2^48, 2^1.51, 2^-1, 1, 2^-45.49, [72, 69, 7+30, 41, 18]),


    'Simon48/96  r17 A': (48, 2^47.5, 2^1.01, 2^-1.5, 1, 2^-45.49, [96, 94, 16+30, 32, 18]),

    'Simon64/96  r17 A': (64, 2^63.5, 2^1, 2^-1.5, 1, 2^-61.5, [96, 92, 24+24, 40, 40]),

    'Simon64/128  r17 A': (64, 2^63.5, 2^1, 2^-1.5, 1, 2^-61.5, [128, 124, 24+42, 40, 21]),

    'Simon96/96  r33 A': (96, 2^95.5, 2^0.4, 2^-1.5, 1, 2^-94.1, [96, 80, 16+16, 80, 80]),
    'Simon96/96  r33 B': (96, 2^95.7, 2^0.6, 2^-1.3, 1, 2^-94.1, [96, 80, 16+16, 80, 80]),

    'Simon96/144  r33 A': (96, 2^95.5, 2^0.4, 2^-1.5, 1, 2^-94.1, [144, 124, 16+30, 80, 66]),

    'Simon128/128  r42 A': (128, 2^127.5, 2^1.52, 2^-1.5, 1, 2^-124.98, [128, 105, 21+25, 107, 103]),

    'Simon128/192  r42 A': (128, 2^127.5, 2^1.52, 2^-1.5, 1, 2^-124.98, [192, 154, 38+25, 90, 103]),

    'Simon128/256  r42 A': (128, 2^127.5, 2^1.52, 2^-1.5, 1, 2^-124.98, [256, 95+113, 38+43, 90, 85]),
}

for d in distinguishers.keys():
    print ("## "+d)
    n, N, lmd_r, lmd_w, a, M, parTime = distinguishers[d]
    print("Parameters: n={} N=2^{:.1f} a={} p_r={:.2f} lambda_R={:.2f} lambda_W={:.2f}\n".format(n, log2(N), a, log2(M), log2(lmd_r), log2(lmd_w) ) )
    F_r = stats.poisson.cdf(a, lmd_r)
    p = 1-F_r
    print ("Succes probability (direct)           : {:.2f}".format(p))


    
    if len(parTime) != 0:

        keysize, K_g, condition, l_0, l_r = parTime

        pairs = log2(N) - 1 + ( n - l_0 )  - l_r

        search = pairs + K_g - condition


        
        print("Parameters: pairs remain={:.2f} keysize K={} Total Guessed Key Bits={}  lambda_W=2^{:.2f}  searchTime {}\n".format( pairs, keysize, K_g, np.log2(lmd_w), search ) )
        F_w = stats.poisson.cdf(a, lmd_w)

        time_es = 2^keysize * (1 - F_w)

        totalTime = N + time_es + 2^search

        #AverageTime = totalTime * (1/p)

        print ("search complexity is (direct)           : 2^{:.2f}\n".format( search ))

        print ("ES time complexity is (direct)           : 2^{:.2f}\n".format(np.log2(time_es) ))

        print ("Total time complexity is          : 2^{:.2f}\n".format(np.log2(totalTime) ))

        #print ("Average time complexity is          : 2^{:.2f}\n".format(np.log2(AverageTime) ))



