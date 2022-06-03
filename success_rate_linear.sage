#! /usr/bin/env sage

# Accompanying code to the paper "Improved Differential and Linear Cryptanalysis on Round-Reduced SIMON"
# Computes the success probabilty of linear attacks with various parameters

log2 = lambda x: log(1.0*x)/log(2.0)
Phi = lambda x: (1 +erf(x/sqrt(2 )))/2 
Phi_inv = lambda x: sqrt(2 )*erfinv(2 *x-1 )


distinguishers = {
    #'Simeck64 29r (5 approx)': (64 , 2^63 , 52 , [2^-58.47 ]+([2^-60.36]*4 )),
    #'Simeck64 29r': (64, 2^62,   26, [2^-58.47]),
    #'Simeck64 30r A': (64, 2^63.5, 24 , [2^-60.36 ]),
    #'Simeck64 30r B': (64, 2^64,   29 , [2^-60.36 ]),
    #'Simeck48 21r': (48, 2^47,   26, [2^-43.56]),
    #'Simeck32 13r': (32, 2^31.5,   37, [2^-27.68]),

    #'Simon128 41r': (128, 2^126, 10, [2^-123.07]),
    #'Simon128 42r': (128, 2^127, 10, [2^-125.07]),
    #'Simon96  33r A': (96, 2^95, 10, [2^-92.60]),
    #'Simon96  33r B': (96, 2^94, 10, [2^-92.60]),



    'Asiacrypt Simon96/144  r33 A': (96, 2^95, 5, [2^-92.61], [144, 130]),


    'Simon48/72  r17 A': (48, 2^47, 10, [2^-45.19], [72, 54]),
    'Simon48/96  r17 B': (48, 2^47, 10, [2^-45.19], [96, 72]),

    'Simon64/96  r23 A': (64, 2^63, 10, [2^-60.25], [96, 79]),

    'Simon64/128  r23 A': (64, 2^63, 10, [2^-60.25], [128, 106]),

    'Simon96/96  r34 A': (96, 2^95, 10, [2^-93.74], [96, 84]),

    'Asiacrypt Simon96/96  r33 B': (10, 2^94, 10, [2^-92.61], [96, 84]),

    'Simon96/144  r33 A': (96, 2^95, 10, [2^-91.74], [144, 130]),

    'Simon96/144  r33 B': (96, 2^94, 10, [2^-91.74], [144, 130]),

    'Simon128/128  r43 A': (128, 2^127, 10, [2^-124.59], [128, 116]),

    'Simon128/256  r43 A': (128, 2^127, 10, [2^-124.59], [256, 230]),

}

for d in distinguishers.keys():

    print ("## "+d)
    n, N, a, M, TimePara = distinguishers[d]

    if n==32:
        # ELP for n==32 is exact
        expC = sum(m for m in M)
        varC = 2*sum(m^2 for m in M)

    else:
        # add 2^-n to ELP
        expC = sum(m+2^-n for m in M)
        varC = 2*sum((m+2^-n)^2 for m in M)
    
    print ("Parameters: n={} N=2^{:.1f} a={} M={} C={:.2f} C'={:.2f}".format(n, log2(N), a, len(M), log2(sum(M)), log2(expC)))

    B = (2^n-N)/(2^n-1)     # Distinct plaintexts
    # B = 1                 # Random plaintext
    
    if len(M) == 1:
        #Single approximation

        sr = B/N+expC
        sw = B/N+2^-n
        
        print(log2(M[0]))
        sN = M[0] / (2^-n) 
        mu = M[0] * N

        p = numerical_approx(2  - 2 *Phi(sqrt(sw/sr)*Phi_inv(1 -2^(-a-1 ))), prec=112 )
        
        #p = numerical_approx(  Phi(  ( sqrt(mu * sN) - Phi_inv(1 - 2^(-a) ) ) / sqrt( sN + 1 )   )   )

        print ("Succes probability (direct)           : {:.2f}".format(p))

        if len(TimePara):
            kG = TimePara[0]
            Time_FWT = TimePara[1]
            attempts = 1/p

            averageTime = (2^(kG-a) + 2^Time_FWT) * attempts

            success_proba = 1 - ( 1 - p )^(attempts) 

            print ("One time complexity (direct)           : {:.2f}".format(log2( (2^(kG-a) + 2^Time_FWT) )))

            print ("Average time complexity (direct)           : {:.2f}".format(log2(averageTime)))

            print ("After {}  attempts, Succes probability rotation (direct)           : {:.2f}".format(attempts, success_proba ) )




    else:
        # Multiple approximations
        mr = B*len(M) + N*expC
        sr = 2*B^2*len(M) + 4*B*N*expC + N^2*varC
        Right = lambda x: RealDistribution('gaussian',sqrt(sr)).cum_distribution_function(x-mr)
    
        w = B+N*2^-n
        Wrong = lambda x: RealDistribution('chisquared',len(M)).cum_distribution_function_inv(x)*w
    
        p = 1 - Right(Wrong(1 - 2^-a))
        print ("Succes probability (gaussian/chi_2)   : {}".format(p))


