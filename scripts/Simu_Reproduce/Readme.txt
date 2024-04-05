To reproduce the results in Results-Sim22_5_1_context, please run 
- simu.R  to store results of Table 1 
- simu_p to store results of Figure 1
- simu_nout to store results of Figure 2
- simu_alpha to store results of Figure 3
- simu_spa to store results of Figure 4
- simu_snr to store results of Figure 5

Other files:
- setup_secure.R for the data generating process
- Function RS1.R for the proposed method
- Function r4.1.R for R^4 method from She and Chen (2017): we slightly modify it to fix the maximum rank of this method for fair comparsion with competing methods
- Function secure.R for sequential cosparse factor regression proposed by Mishra et al. (2017)
- Function RCGL.R for rank constrained group Lasss proposed by Bunea et al. (2012)
- Function RRR.R for reduced rank regression (RRR) for low-rank estimates