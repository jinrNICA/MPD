Quick howto:

1. To write correct DCA into cbmsim tree, use restore_dca
2. Get fitted dca values:
	-- Use get_dca.C to obtain them from cbmsim tree
	-- Use get_fit.C for the 1st iteration of the fit
	-- Use MakeFitDCA.C for final fitted dca values (improve pt efficiency)
3. To get centrality values from multiplicity in TPC use get_centrality
4. To get light data format for the further analysis use create_reduced_tree
5. To get both resolution correction factor and flow use real_flow
