### Figure 3

 - To plot Figure 3, run heatmap.m. To reduce computational expense we executed heatmap.m, in 4 batches: 10 simulations twice and 20 simulations twice. 
 - The results of these simulations are saved in files named batchX_nsim.mat, where X: batch number, and nsim: number of simulations.
 - We processed the batch files using the batch_concat.m script, which aggregated the final parameter values and the percentage of time in the up-state. These processed data are stored in a file named allY.mat, where Y: parameter name.
 - To finally plot the heatmaps, execute the jupyter file Figure 3.ipynb
