import numpy as np
import scipy
import math
from mpmath import *


def compute_def_integral(n_nodrug, py_nodrug, inf_nodrug, stepwidth = 1e-5):
    """compute the distribution of infection incidence
        input: initial population number, infection number, total person-years
        return: the PDF of the infection incidence, the grid of possible rinf values on which the grid was computed....
        """

    mp.dps = 100
    def fun(x):     # function to integrate
        return (power(x, inf_nodrug)) * (power(1 - py_nodrug*x/n_nodrug, n_nodrug-inf_nodrug))
    
    max_r = n_nodrug / py_nodrug
    steps = math.floor(max_r / stepwidth)

    #r_values = np.array([float(i*stepwidth) for i in range(steps)], dtype=float)

    # compute pdf on that grid
    res = []
    for i in range(steps):
        a = i * stepwidth
        b = (i + 1) * stepwidth
        fa = fun(a)
        fb = fun(b)
        res.append((fa + fb) * stepwidth / 2)

    pdf = np.array(res, dtype=float)
    pdf /= pdf.sum()
    return pdf

def compute_efficacy_distribution_bayesian_grid(
    inf_drug, n_drug,
    n_nodrug, py_nodrug,
    n_inf_nodrug, inf_num_sim = None, rinf_null=None, rinf_pdf_sim=None,
    prior=1,
    N_grid=10000,
    flag = 'topdown'
):
    """
    Bayesian efficacy estimation:
    Computes posterior P(E | data_drug, data_control)

    When N_grid == 0, then uses simulation data
    """
    if flag == 'topdown':
        """
        new approach to compute the efficacy distribution using Bayesian inference
        """

        rinf_pdf = compute_def_integral(n_nodrug, py_nodrug, n_inf_nodrug)
        res = np.zeros(101)
        
        r_values = (np.arange(10000)/ 10000) / (n_nodrug / py_nodrug) 
        
        for i in range(101): 
            efficacy = 1 - i / 100
            prob = scipy.stats.binom.pmf(inf_drug, n_drug, r_values * efficacy) 
            res[i] = np.dot(prob, rinf_pdf) 
        

    else: 

        efficacy = 1 - np.arange(101) / 100  # shape (101,)

        # Create a matrix of shape (M, 101) for rinf_null * efficacy
        
        prob_matrix = np.array([
            scipy.stats.binom.pmf(inf_num_sim[j], n_drug, rinf_null[j] * efficacy)
            for j in range(len(rinf_null))
        ])  

        
        res = prob_matrix.sum(axis=0)

    res = np.array(res) * prior

    return res / res.sum()


def sample_incidence(undet_data, trial_name, sample_n):
    # Get trial undet data for the incidence
    trial_undet = undet_data[trial_name]
    n_tot_inf = trial_undet['individuals'] #number of individuals in the study arm
    py_inf = trial_undet['avg_obs_time']*n_tot_inf #total observation time
    n_inf_inf = trial_undet['n_inf'] #number of observed infections
    #Compute pdf for the infection rate)
    pdf_rinf= compute_def_integral(n_tot_inf, py_inf, n_inf_inf)
    # Convert mpmath.mpf values to float64
    pdf_rinf = np.array([float(val) for val in pdf_rinf], dtype=np.float64)
    #Define the interval for your infection rate
    max_rinf = n_tot_inf / py_inf
    step_width = 1e-5
    incidence_values = np.arange(0, max_rinf, step_width)
    # Sample from the incidence values using the probability distribution
    sampled_incidence = np.random.choice(incidence_values[:len(pdf_rinf)], size=sample_n, p=pdf_rinf)
    indices = ((sampled_incidence/ step_width).astype(int))
    pdf_at_samples = pdf_rinf[indices]
    return sampled_incidence, pdf_at_samples

def simulation_individual_infections(n_individuals, py_followup, r_infection_incidence, phi, nsim, trialname, hypo):
    """
    Simulate infection events for each individual with individual efficacy phi.
    Each simulation uses a different infection incidence rate.
    Parameters:
        n_individuals: int
        py_followup: float, average follow-up time (years)
        r_infection_incidence: ndarray (nsim,), baseline infection rate per simulation (/year)
        phi: ndarray (n_individuals,), individual efficacy values (fixed across sims)
        nsim: int, number of simulations

    Returns:
        n_infections: ndarray (nsim,), number infected per simulation
        sim_incidence: ndarray (nsim,), simulated incidence rates
        ref_incidence: ndarray (nsim,), reference incidence rates
    """
    phi = np.asarray(phi)
    r_infection_incidence = np.asarray(r_infection_incidence)
    # print('Simulation shape info:', phi.shape, r_infection_incidence.shape)
    assert phi.shape == (nsim, n_individuals), "phi must be shape (n_individuals,)"
    assert r_infection_incidence.shape == (nsim,), "r_infection_incidence must be shape (nsim,)"
    phi_matrix = phi
    sim_tot_t = 0
    # Compute individual infection rates per simulation
    r_inf_matrix = r_infection_incidence[:, np.newaxis] * (1 - phi_matrix)  # (nsim, n_individuals)
    r_dropoff = 1 / py_followup - r_inf_matrix

    if np.any(r_dropoff < 0):
        raise ValueError("Negative dropout rate detected. Adjust py_followup or r_infection_incidence.")

    # Total rate for competing events (infection or dropout)
    total_rate = r_inf_matrix + r_dropoff  # (nsim, n_individuals)

    # Simulate time to event (infection or dropout)
    tau = np.random.exponential(1 / total_rate)  # (nsim, n_individuals)
    # Simulate event type: infection vs dropout
    rand_event = np.random.random(size=(nsim, n_individuals))

    infected = (tau <= py_followup) & (rand_event < (r_inf_matrix / total_rate))
    n_infections = np.sum(infected, axis=1)  # infections per simulation

    return n_infections, r_inf_matrix


def run_sim_sampled_phi(det_data, undet_data, trial_name, phi_data, p_data, nsim, incidence, pdf_at_samples, hypo):
    trial_det= det_data[trial_name] #get trial data
    inf_num = trial_det['n_inf']
    n_tot = trial_det['individuals'] #number of individuals in the study arm

    trial_undet= undet_data[trial_name] #get trial data
    inf_num_undet = trial_undet['n_inf']
    n_tot_undet = trial_undet['individuals'] #number of individuals in the study arm
    py = trial_undet['avg_obs_time']*n_tot #total observation time

    phi_samples = np.empty((nsim, n_tot))

    for i in range(nsim):
        phi_samples[i, :] = sample_phi_mean(p_data, phi_data, n_tot)

    py_avg = trial_det['avg_obs_time'] 

    infection_reshape, inc_sim = simulation_individual_infections(n_tot, py_avg, incidence, phi_samples, nsim, trial_name, hypo)

    estimated_efficacy = compute_efficacy_distribution_bayesian_grid(inf_num, n_tot,
                             n_tot_undet, py, inf_num_undet, infection_reshape, incidence, pdf_at_samples, flag='bottomup')

    sample_phi = np.random.choice([i for i in range(101)], size=int(1e4), p=estimated_efficacy) #basically sample from the posterior
    
    return infection_reshape, sample_phi


def sample_phi_mean(p, phi, N=100000):
    '''
    p: probability of detection (IPERGAY, iPrEx) or probability adherence (HPTN 083, PURPOSE 2)
    phi: (7, 1000) array with mean efficacy for each individual and adherence strata
    N: number of efficacies to sample, following n_individuals in the trial to simulate
    '''
    num_strata, num_cols = phi.shape #(7, 1000)

    # Sample rows according to p
    row_choice = np.random.choice(num_strata, size=N, p=p)
    col_choice = np.random.randint(0, num_cols, size=N)  # uniform over IDs

    sampled_phi = phi[row_choice, col_choice]

    return sampled_phi

def sample_phi_topdown(p_phi, N):
    '''
    p: probability of efficacy i, from topdown new method
    phi: (101,) array with p of efficacy 0 to 100 
    N: number of efficacies to sample, following n_individuals in the trial to simulate
    '''
    sampled_phi = np.random.choice([i for i in range(101)], size=int(N), p=p_phi)
    
    return sampled_phi


def run_sim_topdown(trial, det_data, undet_data, phi_pdf, n_sim):
    '''
    num_phi_samples: N_individuals
    sample_phi: sample randomly from the distributions phi in phi_data (one per trial)
    
    Important to remember: for each sampled efficacy, 
    N sampled rinf will also be sampled and scaled with the sampled efficacy

    num_phi_samples = 500
    sample_phi = np.random.choice([i for i in range(101)], size=int(num_phi_samples), p=phi)
    num_simulations_to_run = 500 # Number of simulations to run for each sampled efficacy
    '''
    
    det_arm= det_data[trial]
    n_tot = det_arm['individuals'] #number of individuals in the study arm
    inf_num = det_arm['n_inf']
    

    
    # Sample from the incidence values using the probability distribution
    incidence_N = n_sim
    
    sampled_incidence, _ = sample_incidence(undet_data, trial, incidence_N)
    
    phi_samples = np.empty((n_sim, n_tot))
    for i in range(n_sim):
            phi_samples[i, :] = sample_phi_topdown(phi_pdf, n_tot)/100
            #print()
    
    py_avg = det_arm['avg_obs_time'] 
    hypo = 'topdown'
    infection, inc_sim = simulation_individual_infections(n_tot, py_avg, sampled_incidence, phi_samples, n_sim, trial,hypo)
    return infection

def run_hypotheses(
    trialname,
    hypotheses,
    nsim,
    trial_names,
    hyp_list,
    p_sampling,
    phi_topdown,
    phi_,
    trial_parameters_drug_det_men,
    trial_parameters_drug_undet_men,
    seed=123,
):
    """
    Run top-down vs bottom-up for one or multiple hypotheses.

    hypotheses:
        - str: single hypothesis name
        - list[str]: multiple hypotheses
        - "all": run all hypotheses
    """

    np.random.seed(seed)

    # --- normalize hypotheses input ---
    if hypotheses == "all":
        hyp_to_run = hyp_list
    elif isinstance(hypotheses, str):
        hyp_to_run = [hypotheses]
    else:
        hyp_to_run = list(hypotheses)

    # --- trial index ---
    j = trial_names.index(trialname)
    p_det = p_sampling[j]
    phi_tpdwn = phi_topdown[j]
    print(f'\nSimulating {trialname}\n')
    # --- TOP-DOWN once ---
    _, ninf_td = output_topdown(
        trialname,
        trial_parameters_drug_det_men,
        trial_parameters_drug_undet_men,
        phi_tpdwn,
        nsim
    )

    results = {}

    # --- Hypothesis-driven  ---
    print(f'\nHypothesis-driven simulations')
    for hyp in hyp_to_run:
 
        phi_m = phi_[hyp]

        print(f"hypothesis: {hyp}")

        _, ninf_bu, efficacy = output_bottomup(
            trialname,
            trial_parameters_drug_det_men,
            trial_parameters_drug_undet_men,
            p_det,
            phi_m,
            nsim,
            hyp
        )

        pval = p_value(ninf_bu, ninf_td, nsim)

        results[hyp] = {
            "ninf_topdown": ninf_td,
            "ninf_bottomup": ninf_bu,
            "efficacy_posterior": efficacy,
            "p_value": pval,
        }

    return results

def output_topdown(name, trial_det, trial_undet, phi_TD, N):      
    # Topdown part
    print(f'Data-driven simulations')
    infections_tpdwn = run_sim_topdown(name, trial_det, trial_undet, phi_TD, N)
    median, q_025, q_975 = round(np.mean(infections_tpdwn)),np.percentile(infections_tpdwn, 2.5),np.percentile(infections_tpdwn, 97.5)
    print(f'NInf (95% CI): {median} ({round(q_025)}, {round(q_975)})')

    psamp = np.random.choice([i for i in range(101)], size=int(1e4), p=phi_TD)
    eff_mean = round(np.mean(psamp),1)
    eff_low = round(np.percentile(psamp, 2.5),1)
    eff_high = round(np.percentile(psamp, 97.5),1)
    print(f'Efficacy (95% CI): {eff_mean} ({round(eff_low)}, {round(eff_high)})')
    return psamp, infections_tpdwn

def output_bottomup(name, trial_det,trial_undet, p, phi_m, N, hypo):
    incidence_N = N

    incidence, pdf_samples= sample_incidence(trial_undet, name, incidence_N)   
    results = {}

    infections, phi = run_sim_sampled_phi(trial_det, trial_undet, name, phi_m, p, N, incidence, pdf_samples,hypo)

    # infection distribution stats
    ninf_mean = np.mean(infections)
    ninf_low = np.percentile(infections, 2.5)
    ninf_high = np.percentile(infections, 97.5)

    # efficacy posterior stats
    eff_mean = np.mean(phi)
    eff_low = np.percentile(phi, 2.5)
    eff_high = np.percentile(phi, 97.5)
            
    # store results
    if name not in results:
        results[name] = {}
    print('NInf (95% CI):',round(ninf_mean),(round(ninf_low), round(ninf_high)))
    print('Efficacy (95% CI): ', round(eff_mean,1), (round(eff_low,1), round(eff_high,1)))

    results[name][hypo] = {
        "ninf_mean": round(ninf_mean),
        "ninf_ci_low": round(ninf_low),
        "ninf_ci_high": round(ninf_high),
        "eff_mean": round(eff_mean, 1),
        "eff_ci_low": round(eff_low, 1),
        "eff_ci_high": round(eff_high, 1),
    }
    return results, infections, phi

def p_value(infections_bottomup, infections_topdown, sim_num):
    pvalue = np.sum(np.array(infections_bottomup)<=np.array(infections_topdown))/sim_num
    print(f'P-value: {pvalue:.15f}\n')
    return pvalue
