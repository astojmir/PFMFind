"""
This module contains the functions to generate a random 
sequence based on letter composition.

Classes:
    - RandFrag  
        - This class contains the functions which will 
	generate a random sequence of letters based on 
	letter composition.
	

"""

from random import seed, random
from bisect import bisect
from string import join

# Simple random sequence generation
# based on letter composition

class RandFrag:
    """
    This class contains the functions which will generate a random
    sequence of letters based on letter composition.

    """

    def __init__(self, probs_dict):
        """
        Constructor, creates dictionary of probabilities. Initializes 
        a list containing all of the dictionaries (B{cum_probs}).
        Sets the last item in the list to be exactly 1.0.
    
        @param probs_dict: A dictionary of probabilities.  
        """
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
        """
        Generate the random sequence of letters.
    
        @param length: Length of random sequence.
        """
        slst = [self.alphabet[bisect(self.cum_probs, random())] for k in range(length)]
        return join(slst, "")
        
        
        
