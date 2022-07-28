import numpy as np
import argparse
import pandas as pd
import multiprocessing as mp
import time
import random
from numba import jit

@jit(nopython=True)
def set_grid(num_sus, num_s1, num_s2, num_agents):
    if (num_sus+num_s1+num_s2) <= num_agents and num_agents**.5 % 1 == 0:
        frac_sus = num_sus/num_agents
        frac_s1 = num_s1/num_agents
        frac_s2 = num_s2/num_agents
        frac_empty = 1 - frac_sus - frac_s1 - frac_s2
        grid = np.random.choice(np.array([0,1,2,3]), num_agents, [frac_empty,frac_sus,frac_s1,frac_s2])
        return grid
    else:
        print('Error. Number of agents exceeds grid size.')


                                                                  
@jit(nopython=True)
def get_neighbors(index, grid):
    sideLength = int(np.sqrt(len(grid)))
    above = (index - sideLength) % len(grid)
    right = np.floor(index / sideLength)*sideLength + (index % sideLength + 1) % sideLength
    bottom = (index + sideLength) % len(grid)
    left = np.floor(index / sideLength)*sideLength + (index % sideLength - 1) % sideLength
        
    top_right = (right - sideLength) % len(grid)
    bottom_right = (right + sideLength) % len(grid)
    bottom_left = (left + sideLength) % len(grid)
    top_left = (left - sideLength) % len(grid)
    
    res = [above,right,bottom,left,top_right,bottom_right,bottom_left,top_left]
    res = np.array(res).T.astype(np.int64)
    rand_idx = np.random.randint(0, res.size)
    return res[rand_idx]

@jit(nopython=True)
def agent_rules(grid, beta01, beta02, beta12, beta21, v1, v2, prb, eps):
    rand_agent_idx = np.random.randint(0, grid.size)
    agent_state = grid[rand_agent_idx]
    event_id = np.random.randint(1, 9)
    
    if event_id == 1:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < eps:
            grid[rand_agent_idx] = nbr_state
            grid[rand_nbr] = agent_state
    elif event_id == 2 and agent_state == 1:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < prb and nbr_state == 0:
            grid[rand_nbr] = 1
    elif event_id == 3 and agent_state == 2:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < beta01 and nbr_state == 1:
            grid[rand_nbr] = 2
    elif event_id == 4 and agent_state == 3:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < beta02 and nbr_state == 1:
            grid[rand_nbr] = 3
    elif event_id == 5 and agent_state == 2:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < beta21 and nbr_state == 3:
            grid[rand_nbr] = 2
    elif event_id == 6 and agent_state == 3:
        rand_nbr = np.random.randint(0, grid.size)
        nbr_state = grid[rand_nbr]
        if np.random.uniform(0,1) < beta12 and nbr_state == 2:
            grid[rand_nbr] = 3
    elif event_id == 7 and agent_state == 2:
        if np.random.uniform(0,1) < v1:
            grid[rand_agent_idx] = 0
    elif event_id == 8 and agent_state == 3:
        if np.random.uniform(0,1) < v2:
            grid[rand_agent_idx] = 0
            


def sim_run(num_sus, num_s1, num_s2, num_row, num_col, num_steps, beta01, beta02, beta12, beta21, v1, v2, prb, eps, num_sim):
    
    np.random.seed(random.randint(0, 2**32 - 1))

    num_agents = num_row*num_col
    total_steps = num_agents*num_steps

    s_counts = np.zeros(num_steps)
    i1_counts = np.zeros(num_steps)
    i2_counts = np.zeros(num_steps)
    time_steps = np.zeros(num_steps)

    grid = set_grid(num_sus, num_s1, num_s2, num_agents)

    for i in range(total_steps):
        agent_rules(grid, beta01, beta02, beta12, beta21, v1, v2, prb, eps)
        if i % num_agents == 0:
            s_counts[int(i/num_agents)] = np.sum(grid == 1)
            i1_counts[int(i/num_agents)] = np.sum(grid == 2)
            i2_counts[int(i/num_agents)] = np.sum(grid == 3)
            time_steps[int(i/num_agents)] = int(i/num_agents)
            
    d = {'time_steps': time_steps, 
         'sus_counts': s_counts, 
         'strain_1_counts': i1_counts,
         'strain_2_counts': i2_counts}
    df = pd.DataFrame(data=d)
    df.to_csv('output_' + str(num_sim) + '.csv')

def main():
    parser = argparse.ArgumentParser(description='arguments')
    parser.add_argument('num_sus', type=int)
    parser.add_argument('num_s1', type=int)
    parser.add_argument('num_s2', type=int)
    parser.add_argument('num_row', type=int)
    parser.add_argument('num_col', type=int)
    parser.add_argument('num_steps', type=int)
    parser.add_argument('beta01', type=float)
    parser.add_argument('beta02', type=float)
    parser.add_argument('beta21', type=float)
    parser.add_argument('beta12', type=float)
    parser.add_argument('v1', type=float)
    parser.add_argument('v2', type=float)
    parser.add_argument('prb', type=float)
    parser.add_argument('eps', type=float)
    parser.add_argument('n_proc', type=int)
    parser.add_argument('n_sims', type=int)
    args = parser.parse_args()
    
    num_sus = args.num_sus
    num_s1 = args.num_s1
    num_s2 = args.num_s2
    num_row = args.num_row
    num_col = args.num_col
    num_steps = args.num_steps
    beta01 = args.beta01
    beta02 = args.beta02
    beta21 = args.beta21
    beta12 = args.beta12
    v1 = args.v1
    v2 = args.v2
    prb = args.prb
    eps = args.eps
    n_proc = args.n_proc
    n_sims = args.n_sims

    pool = mp.Pool(processes=n_proc)
    simulations = [pool.apply_async(sim_run, args=(num_sus, num_s1, num_s2, num_row, num_col, num_steps, beta01, beta02, beta12, beta21, v1, v2, prb, eps, i)) for i in range(n_sims)]
    output = [sim.get() for sim in simulations]


if __name__=='__main__':
   main()
