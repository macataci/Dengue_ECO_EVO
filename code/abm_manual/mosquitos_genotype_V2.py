import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
import inspect

#Version con genotipos como clase y mosquitos mueren de acuerdo a logística

# Voy a asumir que mosquitos bebés de I o S, ancen infectados. O sea no hay transmisión vertical.

# Infeccion de mosquito:
# - Mosquito suceptible pica a humano infectado.
# - Queda infectado hasta que se muere
# - Asumo que mosquito infectado puede picar infinita cantidad de veces antes de morir.

# Infeccion de humano:
# - Mosquito infectado pica a humano susceptible.


opt = ["A", "C", "T", "G"]        
class Mosquito:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
        
class Human:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
        
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured, mutation_p, K, r):
        
        self.n_humans = n_humans
        self.n_mosquitoes = n_mosquitoes
        self.init_inf_hum = init_inf_hum
        self.init_inf_mos = init_inf_mos
        self.encounter_p = encounter_p
        self.biting_p = biting_p
        self.daysCured = daysCured
        self.gamma = 1/self.daysCured
        self.mutation_p = mutation_p
        self.K = K
        self.r = r
        self.population_hum = []
        self.population_mos = []
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        self.genotypes = pd.DataFrame(columns=["A", "C", "T", "G"])
        self.initialize_pop_hum()
        self.initialize_pop_mos()

    def initialize_pop_mos(self):

        for i in range(self.n_mosquitoes):
            state = "S"
            mosquito = Mosquito(i, state, "")
            self.population_mos.append(mosquito)
        
        sample_infected = random.sample(self.population_mos, self.init_inf_mos)
        
        for mosquito in sample_infected:
            index = self.population_mos.index(mosquito)
            self.population_mos[index].state = "I"
            self.population_mos[index].genotype = random.choice(opt)
            
    def initialize_pop_hum(self):

        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state, "")
            self.population_hum.append(human)
        
        sample_infected = random.sample(self.population_hum, self.init_inf_hum)
        
        for human in sample_infected:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            self.population_hum[index].genotype = random.choice(opt)
            
        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state, self.population_hum[i].genotype)
            
        self.conteos(0)
        self.genotype_counting(0)
        
    def change_state(self):   
 
        # Si ya estaba infectado, muta con base en una probabilidad (humanos nada más)
        # Si no, se infecta con el genotipo de ese mosquito (mosquitos y humanos)
        
        for human in self.population_hum:           
            
            # Human recovery and mosquito infection
            
            # MOSQUITO INFECTION
            
            if human.state == "I":
                for mosquito in self.population_mos:
                    if mosquito.state == "S":                        
                        #TODO dejo multi, transmission o que aca
                        if random.random() > self.encounter_p:                    
                            # Si el encuentro se da, se tiene otra siguiente p a superar, que está relacionada con el biting rate.
                            if random.random() > self.biting_p:
                                mosquito.state ="I"
                                mosquito.genotype = human.genotype

                
                #HUMAN RECOVERY
                
                # Pensar bien si dejarlo así o tipo si pasan tantos días ya se recupera.  
                      
                if random.random() < self.gamma:
                    human.state= "R"
                
                elif random.random () > self.mutation_p:
                    human.genotype = random.choice(opt)
                                
            # Human infection   
            elif human.state == "S":
                
                # Para cada mosquito va a tener una probabilidad de encuentro, si esta es mayor a cierto valor el encuentro se da.
                # Revisar si hacer mejor por grid
                # No voy a revisar a todos pues porque me importan son los infectados
                
                #TODO self.encounter_p que varíe en el tiempo.
                
                for mosquito in self.population_mos:
                    if mosquito.state == "I" and random.random() > self.encounter_p:
                    # Si el encuentro se da, se tiene otra siguiente p a superar, que está relacionada con el biting rate.
                        if random.random() > self.biting_p:
                            human.state ="I"
                            human.genotype = mosquito.genotype
                            # Si se infecta pues se asume que nadie más lo pica
                            break
    
    
    # Puede modificarlo teniendo transmisión vertical.
    # Revisar bien esto porque sería mejor tener en cuenta a cada mosquito
    
    def vital_dynamics(self):
                        
        # dN/dt = rN*(1-(N/k))
        N = len(self.population_mos)
        mosq = round((self.r*N)*(1-(N/self.K)))
        new = mosq - N
        
        for i in range(new):
            id = N + i
            state = "S"
            mosq = Mosquito(id, state, "")
            self.population_mos.append(mosq)
                
                           
    def update(self, t):
        self.change_state()
        for human in self.population_hum: 
            self.data.loc[len(self.data)] = (t,  human.id, human.state, human.genotype)   
        self.vital_dynamics()
    
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
                    
    def genotype_counting (self, t):
        conteo = self.data.where(self.data["Time_step"]==t)["Genotype"].value_counts()
        genotype = conteo.index.tolist()
        if 'A' in genotype:
            A = conteo['A']
        else: 
            A=0
        if 'C' in genotype:
            C = conteo['C']
        else: 
            C=0
        if 'T' in genotype:
            T = conteo['T']
        else: 
            T=0
        if 'G' in genotype:
            G = conteo['G']
        else: 
            G=0  
            
        self.genotypes.loc[t] = (A, C, T, G )   
        
    def run(self, n_steps):
        for t in range(1, n_steps):
            self.update(t)
            self.conteos(t)
            self.genotype_counting(t)
        return self.data, self.counts, self.genotypes

    
# Create DengueABM instance

# (self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured, mutation_p, K, r)

model = SIRmodel(100, 1000, 1, 1, 0.9, 0.9, 14, 0.2, 1500, 1/10)
# Run simulation
df_data, df_conteos, df_genotypes = model.run(100)


time = range(len(df_conteos))

# Plot results

# Dynamics

figure, axis = plt.subplots(1, 2)

axis[0].plot(time, df_conteos["S"].tolist(), label='Susceptible')
axis[0].plot(time, df_conteos["I"].tolist(), label='Infected')
axis[0].plot(time, df_conteos["R"].tolist(), label='Recovered')
axis[0].set_xlabel('Time Step')
axis[0].set_ylabel('Number of Humans')
axis[0].set_title('Humans dynamics')
axis[0].legend()


axis[1].plot(time, df_genotypes["A"].tolist(), label='A')
axis[1].plot(time, df_genotypes["C"].tolist(), label='C')
axis[1].plot(time, df_genotypes["T"].tolist(), label='T')
axis[1].plot(time, df_genotypes["G"].tolist(), label='G')
#axis[1].fill_between(time, df_genotypes["A"].tolist())
#axis[1].fill_between(time,df_genotypes["C"].tolist())
#axis[1].fill_between(time, df_genotypes["T"].tolist())
#axis[1].fill_between(time, df_genotypes["G"].tolist())
axis[1].set_xlabel('Time Step')
axis[1].set_ylabel('Genotype')
axis[1].set_title('Genotype dynamics')
axis[1].legend()

plt.show()

df_data.to_csv("Humans.csv")



#infect = 0

#for i in model.population_mos:
 #   if i.state == "I":
  #      infect+=1





            
