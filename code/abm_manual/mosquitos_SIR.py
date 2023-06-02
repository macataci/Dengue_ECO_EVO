import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
import inspect


# AUN NO METO GENOTIPOS

# Infeccion de mosquito:
# - Mosquito suceptible pica a humano infectado.
# - Queda infectado hasta que se muere
# - Asumo que mosquito infectado puede picar infinita cantidad de veces antes de morir.

# Infeccion de humano:
# - Mosquito infectado pica a humano susceptible.


#TODO matar mosquitos
#TODO assumption de encounter siendo igual h-m, m-h

class Mosquito:
    def __init__(self, id, state): 
        self.id = id
        self.state = state
        
class Human:
    def __init__(self, id, state): 
        self.id = id
        self.state = state
        
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured):
        
        self.n_humans = n_humans
        self.n_mosquitoes = n_mosquitoes
        self.init_inf_hum = init_inf_hum
        self.init_inf_mos = init_inf_mos
        self.encounter_p = encounter_p
        self.biting_p = biting_p
        self.daysCured = daysCured
        self.gamma = 1/self.daysCured
        self.population_hum = []
        self.population_mos = []
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        self.initialize_pop_hum()
        self.initialize_pop_mos()

    def initialize_pop_mos(self):

        for i in range(self.n_mosquitoes):
            state = "S"
            mosquito = Mosquito(i, state)
            self.population_mos.append(mosquito)
        
        sample_infected = random.sample(self.population_mos, self.init_inf_mos)
        
        for mosquito in sample_infected:
            index = self.population_mos.index(mosquito)
            self.population_mos[index].state = "I"
            
    def initialize_pop_hum(self):

        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state)
            self.population_hum.append(human)
        
        sample_infected = random.sample(self.population_hum, self.init_inf_hum)
        
        for human in sample_infected:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            
        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state)
            
        self.conteos(0)
        
    def change_state(self):   
        for human in self.population_hum:           
            # Human recovery and mosquito infection
            if human.state == "I":
                for mosquito in self.population_mos:
                    if mosquito.state == "S":                        
                        #TODO dejo multi, transmission o que aca
                        if random.random() > self.encounter_p:                    
                            # Si el encuentro se da, se tiene otra siguiente p a superar, que está relacionada con el biting rate.
                            if random.random() > self.biting_p:
                                mosquito.state ="I"
                
                # Pensar bien si dejarlo así o tipo si pasan tantos días ya se recupera.        
                if random.random() < self.gamma:
                    human.state= "R"
                                
            # Human infection   
            elif human.state == "S":
                
                # Para cada mosquito va a tener una probabilidad de encuentro, si esta es mayor a cierto valor el encuentro se da.
                # Revisar si hacer mejor por grid
                # No voy a revisar a todos pues porque me importan son los infectados
                for mosquito in self.population_mos:
                    if mosquito.state == "I" and random.random() > self.encounter_p:
                    # Si el encuentro se da, se tiene otra siguiente p a superar, que está relacionada con el biting rate.
                        if random.random() > self.biting_p:
                            human.state ="I"
                            # Si se infecta pues se asume que nadie más lo pica
                            break
                    

                           
    def update(self, t):
        self.change_state()
        for human in self.population_hum: 
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

# (n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured)

model = SIRmodel(30, 500, 3, 10, 0.9, 0.9, 20)
# Run simulation
df_data, df_conteos = model.run(50)

time = range(len(df_conteos))
# Plot results
plt.plot(time, df_conteos["S"].tolist(), label='Susceptible')
plt.plot(time, df_conteos["I"].tolist(), label='Infected')
plt.plot(time, df_conteos["R"].tolist(),label='Recovered')
plt.xlabel('Time Step')
plt.ylabel('Number of Humans')
plt.legend()
plt.show()
infect = 0
for i in model.population_mos:
    if i.state == "I":
        infect+=1

print(infect)




            
