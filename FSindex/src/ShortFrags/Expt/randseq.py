from random import seed, random
from bisect import bisect
from string import join

# Simple random sequence generation
# based on letter composition

class RandFrag:
    def __init__(self, probs_dict):
        self.probs_dict = probs_dict
        self.alphabet = probs_dict.keys()
        self.probs = probs_dict.values()

        self.cum_probs = []
        sum_p = 0 
        for p in self.probs:
            self.cum_probs.append(p + sum_p)
            sum_p += p
        # Make sure the last item is exactly 1.0
        self.cum_probs[-1] = 1.0

        seed()
        
    def rand_frag(self, length):
        slst = [self.alphabet[bisect(self.cum_probs, random())] for k in range(length)]
        return join(slst, "")
        
        
        
