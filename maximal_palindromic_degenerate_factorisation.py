import sys
import itertools
import pysais
import numpy as np

# To install PySAIS module needed for Suffix Arrays and LCP
# follow README.md instructions on https://github.com/AlexeyG/PySAIS

# class RMQ as taken from
# https://gist.github.com/m00nlight/1f226777a49cfc40ed8f
#
# This is most likely not an efficient implementation of RMQ
# Use as a placeholder
#
# When performing queries on a 0-indexed string
# the range considered is [a, b)

class RMQ:
    def __init__(self, n):
        self.sz = 1
        self.inf = (1 << 31) - 1
        while self.sz <= n: self.sz = self.sz << 1
        self.dat = [self.inf] * (2 * self.sz - 1)
    
    def update(self, idx, x):
        idx += self.sz - 1
        self.dat[idx] = x
        while idx > 0:
            idx = (idx - 1) >> 1
            self.dat[idx] = min(self.dat[idx * 2 + 1], self.dat[idx * 2 + 2])
            
    def query(self, a, b):
        return self.query_help(a, b, 0, 0, self.sz)
            
    def query_help(self, a, b, k, l, r):
        if r <= a or b <= l:
            return sys.maxsize
        elif a <= l and r <= b:
            return self.dat[k]
        else:
            return min(self.query_help(a, b, 2 * k + 1, l, (l + r)>>1),
                        self.query_help(a, b, 2 * k + 2, (l + r) >> 1, r))



# converts string into degenerate sequence
def convert_to_degenerate(sequence, k):
    degenerate_sequence = []

    i = 0
    while i < len(sequence):
        if sequence[i] != "[":
            degenerate_sequence.append([sequence[i]])
        else:
            non_solid = []

            i += 1
            while sequence[i] != "]":
                non_solid.append(sequence[i])
                i += 1
                

            degenerate_sequence.append(non_solid)

        i += 1

    return degenerate_sequence



# factorising function
def get_maximal_palindromic_factorisation(degenerate_sequence, k, alphabet):

    lambda_sequence = ""
    lambda_value = 0
    lambda_locations = []
    alphabet_value = {}

    # assign numerical values to alphabet characters
    for i in range(len(alphabet)):
        alphabet_value[alphabet[i]] = i

    print(alphabet_value)

    # replace non solid symbols with lambda characters
    for i in range(len(degenerate_sequence)):
        char = degenerate_sequence[i]

        if len(char) == 1:
            lambda_sequence = lambda_sequence + char[0]
        else:
            lambda_locations.append(i)
            lambda_sequence = lambda_sequence + str(lambda_value)
            lambda_value += 1



    # MATCH TABLE FORMAT:

    # match_table[i][j]
    # i = l1...lk + s1...
    # j = l1...lk

    ##########################
    #   l1 l2 ... lk s1 ...
    # l1
    # l2
    # ...
    # lk
    # s1
    # ...

    # store match table
    match_table = []
    for i in range(len(lambda_locations) + len(alphabet)):
        match_table.append([False] * len(lambda_locations))

    for i in range(len(lambda_locations)):
        for j in range(len(lambda_locations)):
            p1 = 0
            p2 = 0

            char1 = degenerate_sequence[lambda_locations[i]]
            char2 = degenerate_sequence[lambda_locations[j]]

            while p1 < len(char1) and p2 < len(char2):
                if char1[p1] == char2[p2]:
                    print(i, j)
                    match_table[i][j] = True
                    break
                elif char1[p1] < char2[p2]:
                    p1 += 1
                else:
                    p2 += 1

    for i in range(len(alphabet)):
        for j in range(len(lambda_locations)):
            if alphabet[i] in degenerate_sequence[lambda_locations[j]]:
                match_table[len(lambda_locations) + i][j] = True




    # create match checking funtion
    def real_match(char1, char2):           #CHANGE THIS SO THAT LAMBDAS ARE CONSTANT TIME IDENTIFIED
        if char1 in alphabet:
            if char2 in alphabet:
                return char1 == char2
            else:
                return match_table[len(lambda_locations) + alphabet_value[char1]][int(char2)]
        else:
            if char2 in alphabet:
                return match_table[len(lambda_locations) + alphabet_value[char2]][int(char1)]
            else:
                return match_table[int(char1)][int(char2)]


    # create string for suffix tree
    suffix_tree_string = lambda_sequence + '#' + lambda_sequence[::-1] + '$'

    print("lambda_sequence", lambda_sequence)
    print("suffix_tree_string", suffix_tree_string)
    print("match_table", match_table)

    # store length of suffix tree string
    n_st = len(suffix_tree_string)

    # store length of sequence
    n = len(sequence)

    ########################

    # create suffix array and LCP array
    sa = pysais.sais(suffix_tree_string)
    lcp, lcp_lm, lcp_mr = pysais.lcp(suffix_tree_string, sa)

    # create inverse suffix array
    isa = [0] * n_st
    for i in range(0, n_st):
        isa[sa[i]] = i

    # preprocess LCP array for RMQ
    rmq = RMQ(n_st)
    for i in range(n_st):
        rmq.update(i, lcp[i])

    # optional printing of suffix array
    print_suffix_array = False

    if (print_suffix_array):
        for off in sa :
            print('%3d : %s' % (off, suffix_tree_string[off:]))

    ########################

    # define generic LCP function
    def LCP(i, j):
        if i == j:
            return n_st - i

        a = isa[i]
        b = isa[j]

        if a > b:
            a, b = b, a

        return rmq.query(a, b)

    # CHECK REAL_MATCH WORKING - THEN CONTINUE WITH real_LCP
    # Once maximal palindromes are being found, do the graph setup

    ########################

    """
    def real_LCP(i, j, k):
        real_lcp = 0
        mismatch_count = 0

        # perform up to k + 1 LCP queries for each suffix
        while mismatch_count < k + 1:

            # perform an LCP and update next location to check
            real_lcp = real_lcp + LCP(i + real_lcp, j + real_lcp) + 1

            # stop if end of string is reached
            if real_lcp > n - i or real_lcp > n - j:
                break

            # stop if real mismatch encountered
            char1 = lambda_sequence[i + real_lcp]
            char2 = lambda_sequence[j + real_lcp]

            if not real_match(char1, char2):
                return real_lcp - 1

            # increment mismatch count
            mismatch_count += 1

        # remove 1 to ensure final mismatch is not included
        return real_lcp - 1
    
    even_palindromes = [(0, 0)] * (n - 1)
    odd_palindromes = [(0, 0)] * n

    """

    for i in range(len(lambda_locations)):
        for j in range(len(lambda_locations)):
            print(i, j, real_match(str(i), str(j)))
            print(j, i, real_match(str(j), str(i)))
        for j in range(len(alphabet)):
            print(i, alphabet[j], real_match(str(i), alphabet[j]))
            print(alphabet[j], i, real_match(alphabet[j], str(i)))

    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            print(alphabet[i], alphabet[j], real_match(alphabet[i], alphabet[j]))
            print(alphabet[j], alphabet[i], real_match(alphabet[j], alphabet[i]))

    return None


# perform tests
line_count = 0
correct_result_count = 0

for line in open('tests.txt'):

    values = str.split(line)

    k = int(values[0])
    alphabet = values[1]
    sequence = values[2]
    correct_result = values[3]

    print(k, sequence)

    degenerate_sequence = convert_to_degenerate(sequence, k)
    print(degenerate_sequence)
    factorisation = get_maximal_palindromic_factorisation(degenerate_sequence, k, alphabet)

    # print result of each test
    if (factorisation == None):
        print("NO FACTORISATION")
    else:
        print("FACTORISATION: ", factorisation)

    # increment test result counts
    if factorisation == correct_result:
    	correct_result_count += 1
    	print("Correct")
    else:
    	print("Incorrect")
    print("")

    line_count += 1


print(str(correct_result_count) + " out of " + str(line_count) + " tests passed")
	
