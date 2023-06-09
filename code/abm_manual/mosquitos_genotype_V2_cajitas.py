import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
import inspect
from statistics import mean 
import time

# Lo mismo de V2 pero ahora quiero que los mosquitos no sean agentes sino cajas.
# Ahora el genoma que tiene cada humano es el que más se repite.


# Voy a asumir que mosquitos bebés de I o S, ancen infectados. O sea no hay transmisión vertical.
# TODO revisar esto a futuro
# Infeccion de mosquito:
# - Mosquito suceptible pica a humano infectado.
# - Queda infectado hasta que se muere
# - Asumo que mosquito infectado puede picar infinita cantidad de veces antes de morir.

# Infeccion de humano:
# - Mosquito infectado pica a humano susceptible.


opt = ["A", "C", "T", "G"]        

class Human:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
        
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, hum_t_inf, mutation_p, K, r, mosq_t_inf):
        
        self.n_humans = n_humans
        self.n_mosquitoes = n_mosquitoes
        self.init_inf_hum = init_inf_hum
        self.init_inf_mos = init_inf_mos
        self.encounter_p = encounter_p
        self.biting_p = biting_p
        self.hum_t_inf = hum_t_inf
        self.gamma = 1/self.hum_t_inf
        self.mutation_p = mutation_p
        self.mosq_t_inf = 1/mosq_t_inf
        self.K = K
        self.r = r
        self.population_hum = []
        self.population_mos_S = self.n_mosquitoes-self.init_inf_mos
        self.population_mos_I = self.init_inf_mos
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        self.genotypes = pd.DataFrame(columns=["A", "C", "T", "G"])
        self.initialize_pop_hum()


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
            if human.state == "I": 
                for mosquito in range(self.population_mos_S):    
                    if random.random() > self.encounter_p:                    
                        if random.random() > self.biting_p:
                            self.population_mos_I+=1
                            self.population_mos_S-=1
                                
                if random.random () > self.gamma:
                    human.state = "R"
                
                elif random.random () > self.mutation_p:
                    human.genotype = random.choice(opt)
                                    
        # Human infection   
            elif human.state == "S":
                for mosquito in range(self.population_mos_I):
                    if random.random() > self.encounter_p:
                        if random.random() > self.biting_p:
                            human.state ="I"
                            human.genotype = self.data["Genotype"].value_counts().idxmax()
                            # Si se infecta pues se asume que nadie más lo pica
                            break
                    
        for mosquito in range(self.population_mos_I):
            if random.random() > 1/self.mosq_t_inf:
                self.population_mos_I-=1
                    
    # Puede modificarlo teniendo transmisión vertical.
    # Revisar bien esto porque sería mejor tener en cuenta a cada mosquito
    
    def vital_dynamics(self):
                        
        # dN/dt = rN*(1-(N/k))
        N = self.population_mos_I + self.population_mos_S
        mosq = round((self.r*N)*(1-(N/self.K)))
        self.population_mos_S+=mosq  
                           
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
        #mosquitos = [len(self.population_mos)]
        for t in range(1, n_steps):
            self.update(t)
            self.conteos(t)
            self.genotype_counting(t)
           # mosquitos.append(len(self.population_mos))
        return self.data, self.counts, self.genotypes

start_time = time.time()
# (self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured, mutation_p, K, r)
sims = 2
dias = 60
estados = 3
#x,y,z = simulaciones, tiempos, estado
matriz = np.zeros((sims, dias, estados))
for i in range(sims):
    model = SIRmodel(100, 1000, 1, 200, 0.95, 0.95, 14, 0.2, 1500, 1/10, 4)
    df_data, df_conteos, df_genotypes = model.run(dias)
    S = df_conteos["S"].tolist()
    I = df_conteos["I"].tolist()
    R = df_conteos["R"].tolist()
    matriz[i, :, 0] = S
    matriz[i, :, 1] = I
    matriz[i, :, 2] = R

# Plot results
mediaS_total = []
mediaI_total = []
mediaR_total = []
quantiles = [50, 95]

statistics_50 = pd.DataFrame(columns=["LowS", "HighS", "LowI", "HighI", "LowR", "HighR"])
statistics_95 = pd.DataFrame(columns=["LowS", "HighS", "LowI", "HighI", "LowR", "HighR"])
low_q50 = ((100-50)/2)/100
high_q50 = 1-low_q50

low_q95 = ((100-95)/2)/100
high_q95 = 1-low_q95

for t in range(dias):
    S=matriz[:, t, 0]
    I=matriz[:, t, 1]
    R=matriz[:, t, 2]
    mediaS_t = round(mean(S))
    mediaS_total.append(mediaS_t)
    mediaI_t = round(mean(I))
    mediaI_total.append(mediaI_t)
    mediaR_t = round(mean(R))
    mediaR_total.append(mediaR_t)
    statistics_50.loc[t] = (np.quantile(a= S, q=low_q50), np.quantile(a= S,q=high_q50), np.quantile(a= I, q=low_q50), np.quantile(a= I, q=high_q50), np.quantile(a= R,q=low_q50), np.quantile(a= R,q=high_q50))
    statistics_95.loc[t] = (np.quantile(a= S, q=low_q95), np.quantile(a= S,q=high_q95), np.quantile(a= I, q=low_q95), np.quantile(a= I, q=high_q95), np.quantile(a= R,q=low_q95), np.quantile(a= R,q=high_q95))

# Dynamics
colors = ["#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

figure, axis = plt.subplots(1, 2)
times =range(dias)
axis[0].plot(times, mediaS_total, label="S", color=colors[0])
axis[0].fill_between(times, statistics_95["LowS"], statistics_95["HighS"], alpha=0.2, color=colors[0])
axis[0].fill_between(times, statistics_50["LowS"], statistics_50["HighS"], alpha=0.3, color=colors[0])

axis[0].plot(times, mediaI_total, label="I", color=colors[1])
axis[0].fill_between(times, statistics_95["LowI"], statistics_95["HighI"], alpha=0.2, color=colors[1])
axis[0].fill_between(times, statistics_50["LowI"], statistics_50["HighI"], alpha=0.3, color=colors[1])

axis[0].plot(times, mediaR_total, label="R", color=colors[2])
axis[0].fill_between(times, statistics_95["LowR"], statistics_95["HighR"], alpha=0.2, color=colors[2])
axis[0].fill_between(times, statistics_50["LowR"], statistics_50["HighR"], alpha=0.3, color=colors[2])

axis[0].set_xlabel('Time Step')
axis[0].set_ylabel('Number of Humans')
axis[0].set_title('Humans dynamics')
axis[0].legend()


axis[1].plot(times, df_genotypes["A"].tolist(), label='A')
axis[1].plot(times, df_genotypes["C"].tolist(), label='C')
axis[1].plot(times, df_genotypes["T"].tolist(), label='T')
axis[1].plot(times, df_genotypes["G"].tolist(), label='G')
#axis[1].fill_between(time, df_genotypes["A"].tolist())
#axis[1].fill_between(time,df_genotypes["C"].tolist())
#axis[1].fill_between(time, df_genotypes["T"].tolist())
#axis[1].fill_between(time, df_genotypes["G"].tolist())
axis[1].set_xlabel('Time Step')
axis[1].set_ylabel('Genotype')
axis[1].set_title('Genotype dynamics')
axis[1].legend()

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))






            
