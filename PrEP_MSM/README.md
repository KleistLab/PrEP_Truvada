# PrEP Simulation for MSM

## 1. Requirements

### Software dependencies
```
mpmath==1.3.0
numpy==2.4.0
pandas==2.3.3
scipy==1.16.3
torch==2.2.0
```
### Operating systems
- macOS Ventura 13.5.1

### Hardware
- Cluster or high-performance computing environment required for mechanistic hypothesis-driven simulations involving full cohort of virtual patients

### Instructions

1. Clone the repository:
```
git clone https://github.com/KleistLab/PrEP_Truvada.git
cd PrEP_Truvada/PrEP_MSM/
```
2. Create virtual environment: 
```
python -m venv venv
source venv/bin/activate  # macOS/Linux
```
3. Install packages
Installation time ~5-10 minutes.

```
pip install -r requirements.txt
```
### Demo

1. Hypothesis testing (main results)
- Open the Jupyter notebook notebooks/hypothesis_test_demo.ipynb
- Follow the cells to load trial data and run simulations for the hypothesis outlined.
- Expected simulation time is of ~5 minutes for 100,000 simulations per trial for all hypothesis on a standard desktop.
