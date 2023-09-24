#%%
"To run the simulation"
sim = S.Simulation(10, 1, 0.05, [0,0], [0,0] , con_R = 10)
sim.run(1000, True, False, False, False)