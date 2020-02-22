import random as rn

def HMM_generator(p_o,e_tau,m_tau,pi,n):
    
    # intialize state and output sequences
    state_seqs, output_seqs = "", ""

    # add first states to state sequences and output

    # 1. determine min and max probabilities for checking
    min_p, max_p = min(pi.keys()), max(pi.keys())

    # 2. add state given a random throw
    state_seqs = pi[min_p] if round(rn.random(),2) < min_p else pi[max_p]

    # 3. add output given another random throw
    min_p, max_p = min(p_o.keys()), max(p_o.keys())
    output_seqs = p_o[min_p] if round(rn.random(),2) < min_p else p_o[max_p]

    # 4. decrement n
    n -= 1
    
    # generate the remainder of the sequences

    # iterate through the remaining range
    for i in range(n):

        # Case 1: Last state visited is E
        if state_seqs[len(state_seqs)-1] == 'E':

            # p(e) tau probabilities
            min_tau, max_tau = min(e_tau.keys()),max(e_tau.keys())

            # through dice, go to new state
            current_state = e_tau[min_tau] if round(rn.random(),2) < min_tau else e_tau[max_tau]
            state_seqs += current_state

        # Case 2: Last state visited is M
        else :
            
            # p(m) tau probabilities
            min_tau, max_tau = min(m_tau.keys()),max(m_tau.keys())

            # through dice, go to new state
            current_state = m_tau[min_tau] if round(rn.random(),2) < min_tau else m_tau[max_tau]
            state_seqs += current_state

        # add output state
        min_p, max_p = min(p_o.keys()), max(p_o.keys())
        current_output = p_o[min_p] if round(rn.random(),2) < min_p else p_o[max_p]
        output_seqs += current_output

    return state_seqs,output_seqs

# RANDOM SIMULATOR #

def random_simulator(p_o, e_tau, m_tau, pi, i):
    
    for y in range(1,i+1):

        # 1. Initial probabilities - dice for initial probabilities
        dice = round(rn.random(),2)
        pi[dice] = 'E'
        pi[round(1-dice,2)] = 'M'

        # 2. Output probabilities - dice for initial probabilities
        dice = round(rn.random(),2)
        p_o[dice] = 'T'
        p_o[round(1-dice,2)] = 'H'

        # 3. P(E) Probabilities
        dice = round(rn.random(),2)
        e_tau[dice] = 'E'
        e_tau[round(1-dice,2)] = 'M'

        # 4. P(M) Probabilities
        dice = round(rn.random(),2)
        m_tau[dice] = 'E'
        m_tau[round(1-dice,2)] = 'M'

        # 5. Generate a random length
        n = rn.randint(1,100)
        
        print(f' Iteration: {y}\n')
        print(f' Initial Probabilities: {pi}')
        print(f' Tau State E: {e_tau}')
        print(f' Tau State M: {m_tau}')
        print(f' Output Probabilities: {p_o}')
        print(f' Random Generated Sequence Length: {n}\n')

        results = HMM_generator(p_o,e_tau,m_tau,pi,n)
        states_ratio = {'%E':round(results[0].count('E')/len(results[0]),2),'%M':round(results[0].count('M')/len(results[0]),2)}
        output_ratio = {'%H':round(results[1].count('H')/len(results[1]),2),'%T':round(results[1].count('T')/len(results[1]),2)}

        print(f'State Results: {results[0]}')
        print(f'{states_ratio}\n')
        print(f'Output Results: {results[1]}')
        print(f'{output_ratio}\n')
        
        p_o, e_tau, m_tau, pi = dict(), dict(), dict(), dict()
        
p_o, e_tau, m_tau, pi = dict(), dict(), dict(), dict()
random_simulator(p_o, e_tau, m_tau, pi,20)


