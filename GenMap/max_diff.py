import math, sys

# This function is adapted from the function bwa_cal_maxdiff (bwtaln.c)
# https://github.com/lh3/bwa/blob/master/bwtaln.c
def bwa_cal_maxdiff(l: int, err: float, thres: float) -> int:
    elambda = math.exp(-l * err)
    
    y = 1.0
    x = 1
    sum_val = elambda
    
    for k in range(1, 1000):
        y *= l * err
        x *= k
        sum_val += elambda * y / x
        
        if (1.0 - sum_val) < thres:
            return k
            
    return 2

if __name__ == "__main__":
    kmer_len = int(sys.argv[1])
    error_rate = float(sys.argv[2])
    threshold = float(sys.argv[3])

    print(bwa_cal_maxdiff(kmer_len, error_rate, threshold)) 