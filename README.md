# Eco-evolutionary tradeoffs in the dynamics of prion strain competition
Code for simulations and analysis for the paper "Eco-evolutionary tradeoffs in the dynamics of prion strain competition"

## Contents

- `npm_ode_analysis`:
  - `npm_equilibrium.nb`: Mathematica code for determining the equilibrium of the NPM model
  - `npm_odes.ipynb`: Python notebook that explores the dynamics of the NPM model
- `stochastic_sims`:
  - `infected_cell_invasion`: Code that simulates inter-cellular prion infection among cells that have already been infected
  - `uninfected_cell_invasion`: Code that simulates inter-cellular prion infection of cells that are uninfected
  - `prion_b_evolution`: Simulation of the evolution of the fragmentation rate among a prion population
- `two_strain_model`:
  - `compartmental_ode`:
    - `compartmental_model_stability.nb`: Mathematica code for the stability analysis of the two-strain model
    - `compartmental_ode_plots.ipynb`: Python notebook that explores the dynamics of the two-strain ODE model
  - `agent_based_model`:
    - `agent_sim_two_strains_global.py`: Python code for agent-based simulation of two strain model with a global neighborhood
    - `agent_sim_two_strains_moore.py`: Python code for agent-based simulation of two strain model with a moore neighborhood
    - `agent_sim_visual.ipynb`: Code that visualizes the dynamics of the two-strain agent based model on a grid


  
