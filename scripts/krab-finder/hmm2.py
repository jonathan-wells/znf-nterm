#!/usr/bin/env python3

import numpy as np
from scipy.special import logsumexp
import matplotlib.pyplot as plt


class HMM(object):
    """Base class for HMMs."""

    def __init__(self, nstates, initprobs, transition, emission):
        self.states = np.array([i for i in range(nstates)])
        self.initprobs = initprobs
        self.transition = transition  # Transition probs between hidden states
        self.emission = emission  # Emission probs given hidden states
        self.emission_states = np.array([i for i in range(emission.shape[1])])
        self._check_input()


    
    def _check_input(self, err=1e-5):
        """Checks transition and emission probability matrix format."""
        if self.initprobs.shape != self.states.shape:
            raise ValueError('incorrect initprob format')
        if self.initprobs.shape[0] != self.transition.shape[0]:
            raise ValueError()
        if self.emission.shape[0] != self.transition.shape[0]:
            raise ValueError()
        if (1 - sum(self.initprobs)) > err:
            raise ValueError('initprobs do not sum to 1')
        for i in self.states:
            if abs(1 - sum(self.transition[i])) > err:
                raise ValueError('transition probs do not sum to 1')
        for i in self.states:
            if abs(1 - sum(self.emission[i])) > err:
                raise ValueError('emission probs do not sum to 1')

    def simulate(self, t):
        """Given t total states, simulate states and emitted data."""
        simhidden = [np.random.choice(self.states, p=self.initprobs)]
        for i in range(1, t):
            pred = np.random.choice(self.states, 
                                    p=self.transition[simhidden[i-1], :])
            simhidden.append(pred)
        
        simobs = []
        for i in simhidden:
            pred = np.random.choice(self.emission_states,
                                    p=self.emission[i, :])
            simobs.append(pred)
        return simobs, simhidden

    def viterbi(self, observations):
        """Calculates the maximum likelihood state path given observed data."""
        log_initprobs = np.log(self.initprobs)
        log_transition = np.log(self.transition)
        log_emission = np.log(self.emission)
    
        tracemat = np.zeros((self.states.shape[0], len(observations)), int)
        scoremat = np.zeros((self.states.shape[0], len(observations)))
        
        # Populate matrices
        for i in self.states:
            scoremat[i, 0] = log_initprobs[i] + \
                             log_emission[i, observations[0]]
        def prob(i, j, k, obs):
            return scoremat[k, j-1] + \
                   log_transition[k, i] + \
                   log_emission[i, obs]
        
        for j, obs in enumerate(observations[1:], 1):
            for i in self.states:
                probs = np.array([prob(i, j, k, obs) for k in self.states])
                scoremat[i, j] = probs.max()
                tracemat[i, j] = np.argmax(probs)
         
        # Traceback 
        z = np.argmax(scoremat[:, -1])
        preds = []
        scores = []
        j = len(observations) - 1
        while j >= 0:
            preds.append(z)
            scores.append(scoremat[z, j])
            z = tracemat[z, j]
            j -= 1
        
        return preds[::-1], scores

    def forward(self, observations):
        # TODO: Tricky bit here is dealing with logs of sums. This problem does
        # not exist with Viterbi because we simply take the max prob at each
        # step. For the forward algorithm, we can either take the computational
        # hit or use log interpolation tables to approximate the correct answer.
        
        log_initprobs = np.log(self.initprobs)
        log_transition = np.log(self.transition)
        log_emission = np.log(self.emission)
        
        pathprob = []

        startprob = []
        for i in self.states:
            x = log_emission[i, observations[0]] + log_initprobs[i]
            startprob.append(x)
        pathprob.append(logsumexp(startprob))

        def prob(i, j, k, obs):
            return pathprob[-1] + \
                   log_transition[k, i] + \
                   log_emission[i, obs]

        for j, obs in enumerate(observations[1:], 1):
            probs = []
            for i in self.states:
                for k in self.states:
                    probs.append(prob(i, j, k, obs))
            pathprob.append(logsumexp(probs))
        
        return np.exp(pathprob)
                
    def backward(self, observations):
        return self.forward(observations[::-1])[::-1]

    def posterior_decoding(self, observations):
        p = self.forward(observations)
        q = self.backward(observations)
        # pathp = np.log()
        post = []
        for i, val in enumerate(p):
            post.append((val*q[i])/p[-1])
        return post
                

def test_viterbi(n):
    state_code = ['G', 'I']
    obs_code = ['A', 'C', 'T', 'G']
    nstates = len(state_code)
    initprobs = np.array([0.99, 0.01], float)
    transition = np.array([[0.99, 0.01],
                           [0.05, 0.95]], float)
    emission = np.array([[0.4, 0.1, 0.4, 0.1],
                         [0.2, 0.3, 0.2, 0.3]], float)
    test = HMM(nstates, initprobs, transition, emission)
    so, sh = test.simulate(n)
    preds, scores = test.viterbi(so)
    matches = ''
    for i, val in enumerate(sh):
        if val == preds[i]:
            matches += '|'
        else:
            matches += ' '
    print('observations:', ''.join([obs_code[i] for i in so]))
    print()
    print('states:      ', ''.join([state_code[i] for i in sh]))
    print('             ', matches),
    print('predictions: ', ''.join([state_code[i] for i in preds]))

def test_forward(n):
    state_code = ['G', 'I']
    obs_code = ['A', 'C', 'T', 'G']
    nstates = len(state_code)
    initprobs = np.array([0.9, 0.1], float)
    transition = np.array([[0.999, 0.001],
                           [0.05, 0.95]], float)
    emission = np.array([[0.05, 0.45, 0.05, 0.45],
                         [0.2, 0.3, 0.2, 0.3]], float)
    test = HMM(nstates, initprobs, transition, emission)
    for i in range(1):
        so, sh = test.simulate(100)
        post = test.posterior_decoding(so)
        # print([f'{i:.2e}' for i in post ])
        plt.plot(np.arange(len(post)), post)
        plt.show()

def amino_acid_attributes():
    aadict = {}
    with open('../../data/amino_acid_attributes.txt') as infile:
        for line in infile:
            line = line.strip().split('\t')
            code, attribute = line[2], line[4]
            aadict[code] = attribute
    return aadict

def main():
    seq = """MKVIVIRHMHFDGDDDAGGADDDDDGDTDYDDDDDGYDDFDDDDEDNDTDTDDGVDIDDD\
DTDDENGSSEGINQPLDAMNSDEHDEVIGDVHVCGNCRGEFVVFADFVKHKQNCIKKQVV\
MIYNEGQEQMEDSSRMSEKDQDTASDRSMTPDSRKSVEPNYVSEGKAQQQDEQPGRSSDQ\
SPNPSLPPSNVSLDNLENTPYAVAQFPRESAMVPGHPDVVMIREQLYALQQQHLHQLHMI\
QSIQHQINHLTAQQQRQQIQHQQLQQLQQQNQPKQEPKDSSGSQESPELDEKPLESKTEP\
APSSDTNTATTSSIPSTPIVPIPTIPGMIHPPHILSPAGAAAASAVQQSAIEAMHAQGRL\
LPPLPPPFHLPSQKSPFEMLQQHANRHMPLNPLGIPPPLPYHESLLSGNQKAKPPNVSVF\
DPKFSADDPFFRHKCRFCHKVFGSDSALQIHIRSHTGERPYKCNICGNRFSTKGNLKVHF\
QRHQSRYPHIHMDVRHHPPTPPLMDQPLQPPSIHQIPISTIPNDQARAMYEEDMLPHQPS\
KSTSPRESEGSMADSSRHYQEEEMKQYNSADRNGDDDHDSPSTSQSGSLSPNERRQTTST\
SPSIALSSNVFPVMSAGGIPMSSGYPLPYTGFPRMAFGQNGIMPFESSKSAADVDPVKPE\
TPPSLANLAEAGMRASETSKLQQLVENIEQKITDPNQCVICHRVLSCKSSLQLHYRTHTG\
ERPFKCKICGRSFTTKGNLKTHYAVHRSKGPIRILHDCPICEKRFSNPLVLQQHLRLHAA\
EGMTIPPHLDHPRSEEMHMKDDEDDIEKTVGEDNHMDVDERNDLPTPEEGELPDDTFPEG\
DKPEDGTLPSDVNSPPAENHASESRASDDEKETGQLQSRPFMPSILGLQPVSSMANGDSS\
TAPSPGHSVDTKSMMSEDENRNINMAMYNHDGPESISSSMADEEYDGNYKNEENDTDMSG\
ALDLRPKGTPQENIPPHSRDTGDMPGYYQHDGHYLEDVKPNMTSAGLLMSVPSMQLGSRR\
PMRPNTTCHYCGKVFACTSALNIHYRSHTKERPFRCDVCSKGFSTKGNLKQHMLTHKIRD\
MPMPRYDSPINIPGNSNAKSPVMPKPISQELPLQQTTTTNSQTPPLKRSLSVSPEEANKR\
SQFKHSCSVCGKQFQTNSALEIHMRTHTGEKPFRCHICSRAFTTKGNLKVHLGTHMWNGG\
GRRGRRIAVEPPFIVSPKEGEIFRGDLFGMHPGFAEGAFYHYPPPHLLNGFAAAAANSKA\
NEISVIQQPTSLHLQDSGSVIMRNSLAPSVESREATYSREAERHREMVLSRDEALNHEIA"""
    aadict = amino_acid_attributes()
    propseq = ''.join([aadict[c] for c in seq])
    # print(seq)
    print(propseq)

if __name__ == '__main__':
    main()
