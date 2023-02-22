#!/usr/bin/env python3

import numpy as np
from scipy.special import logsumexp

class MarkovChain():
    def __init__(self,
                 states,
                 initial_probabilities,
                 transition_probabilities,
                 ) -> None:
        self.state_labels = {label: state for (state, label) in enumerate(states)} 
        self.states = np.array([i for i in range(len(states))], dtype='int')
        self.initial_probabilities = np.array(initial_probabilities, dtype='float32')
        self.transition_probabilities = np.array(transition_probabilities, dtype='float32')
        self._check_input()
        self.initial_probabilities = np.log(self.initial_probabilities)
        self.transition_probabilities = np.log(self.transition_probabilities)

    def _check_input(self, err=1e-2):
        """Checks transition_probabilities and emission probability matrix format."""
        if self.initial_probabilities.shape != self.states.shape:
            raise ValueError('incorrect initprob format')
        if self.initial_probabilities.shape[0] != self.transition_probabilities.shape[0]:
            raise ValueError()
        if (1 - sum(self.initial_probabilities)) > err:
            raise ValueError('initial_probabilities do not sum to 1')
        for i in self.states:
            if abs(1 - sum(self.transition_probabilities[i])) > err:
                print(sum(self.transition_probabilities[i]))
                raise ValueError('transition_probabilities probs do not sum to 1')
    
    def log_likelihood(self, sequence):
        sequence = [self.state_labels[c] for c in sequence]
        log_likelihood = self.initial_probabilities[sequence[0]]
        for i in range(1, len(sequence)):
            log_likelihood += self.transition_probabilities[sequence[i-1], sequence[i]]
        return log_likelihood



def load_PAM30():
    with open('../../data/pam30.txt') as infile:
        infile.readline()
        amino_acids = infile.readline().strip().split('\t')
        transition_probabilities = []
        for line in infile:
            line = [2**(float(i)) for i in line.strip().split()]
            transition_probabilities.append(line)
    return amino_acids, transition_probabilities
    
def main():
    states, transition_probabilities = load_PAM30()
    initial_probabilities = [1.0/len(states)]*len(states)
    # initial_probabilities = [0.25, 0.25, 0.25, 0.25]
    # transition_probabilities = [[0.2, 0.1, 0.4, 0.3],
    #                             [0.1, 0.2, 0.4, 0.3],
    #                             [0.4, 0.1, 0.2, 0.3],
    #                             [0.3, 0.2, 0.3, 0.2]]
    mc = MarkovChain(states, initial_probabilities, transition_probabilities)

    observed_sequence = 'MAFIKEESED'
    ll = mc.log_likelihood(observed_sequence)
    print(ll)

if __name__ == '__main__':
    main()
