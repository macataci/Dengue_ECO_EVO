import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display

#Version hiper sencilla, sin mosquitos y sin contactos entre humanos.
        
class Human:
    def __init__(self, id, state): 
        self.id = id
        self.state = state
        
class SIRmodel:
    def __init__(self, n_humans, init_inf_hum, beta, daysCured):
        self.n_humans = n_humans
        self.init_inf_hum = init_inf_hum
        self.daysCured = daysCured
        self.beta = beta
        self.gamma = 1/self.daysCured
        self.population = []
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        self.initialize_pop()

        
    def initialize_pop(self):

        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state)
            self.population.append(human)
        
        sample_infected = random.sample(self.population, self.init_inf_hum)
        
        for human in sample_infected:
            index = self.population.index(human)
            self.population[index].state = "I"
            
        for i in range(len(self.population)):
            self.data.loc[i]=(0, self.population[i].id, self.population[i].state)
            
        self.conteos(0)
        
    def change_state(self):
        
        #me falta meter a los mosquitos acá je
        #revisar esto bien pr qué podría decir como: quedas infectado 14 días, o solo manejarlo con p de recover.
        for human in self.population:
            
            #Infection
            
            if human.state=="S":
                if random.random()< self.beta:
                    human.state="I"
                    
            #Recovery
            
            elif human.state=="I":
                if random.random()< self.gamma:
                    human.state="R" 
                           
    def update(self, t):
        self.change_state()
        for human in self.population: 
            self.data.loc[len(self.data)] = (t,  human.id, human.state)
    
    def conteos(self, t):
        conteo = self.data.where(self.data["Time_step"]==t)["State"].value_counts()
        state = conteo.index.tolist()
        if 'S' in state:
            S = conteo['S']
        else: 
            S=0
        if 'I' in state:
            I = conteo['I']
        else: 
            I=0
        if 'R' in state:
            R = conteo['R']
        else: 
            R=0
            
        self.counts.loc[t] = (S, I, R)          
                    
    def run(self, n_steps):
        for t in range(1, n_steps):
            self.update(t)
            self.conteos(t)
        return self.data, self.counts

    
# Create DengueABM instance
model = SIRmodel(500, 30, 0.05, 14)

# Run simulation
df_data, df_conteos = model.run(1000)

time = range(len(df_conteos))
# Plot results
plt.plot(time, df_conteos["S"].tolist(), label='Susceptible')
plt.plot(time, df_conteos["I"].tolist(), label='Infected')
plt.plot(time, df_conteos["R"].tolist(),label='Recovered')
plt.xlabel('Time Step')
plt.ylabel('Number of Humans')
plt.legend()
plt.show()
#500,1000-3:46-4:45


            
