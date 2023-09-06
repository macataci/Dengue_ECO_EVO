import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
import inspect
from statistics import mean 
import time

"""
En este script se presentarán las siguientes características:
- Se tendrá un modelo SIR en el cual los humanos se infectan al contacto con mosquitos infectados.
- Igualmente, los mosquitos se infectan al picar a humanos infectados.
- Los mosquitos serán cajas y no agentes dado a que no almacenaremos información relevante de ellos, tipo el genoma.
- El genoma asignado a los humanos saldrá de un pool creado inicialmente. Este pool corresponde a genotipos de un mismo serotipo, por lo que 
    sus distancias genéticas son pequeñas. El genoma asignado será el de mayor frecuencia.
- Se tendrá mutación (pero pequeñita) y cuando esto suceda, esta nueva secuencia se añadirá al pool y se quitará la que se tenía antes de la mutación.
- También en caso de recuperación, ese genoma sale del pool, y se recuperan con una probabilidad uniforme. 
"""


#TODO revisar dinámicas vitales de humanos y mosquitos
#TODO mirar si crear una clase de tipo genotype o como creo ese pool
opt = ["A", "C", "T", "G"]      

class Human:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
               
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, hum_t_inf, mutation_p, K, r, mosq_t_inf):
        
        """
            Args:
                n_humans      (int): human population
                n_mosquitoes  (int): mosquito population
                init_inf_hum  (int): initial number of infected humans          
                init_inf_mos  (int): initial number of infected mosquitoes
                encounter_p (float): probability of encounter of a human with a mosquito and viceverse
                bitting_p   (float): probability that the encounter results in a bite from the mosquito
                hum_t_inf   (float): days that a human is infected
                mutation_p  (float): mutation rate of the genome
                K             (int): carrying capacity mosquitoes
                r           (float): constant of proportionality
                mosq_t_inf  (float): days that a mosquito is infected   
            """
        self.n_humans = n_humans
        self.n_mosquitoes = n_mosquitoes
        self.init_inf_hum = init_inf_hum
        self.init_inf_mos = init_inf_mos
        self.encounter_p = encounter_p
        self.biting_p = biting_p
        self.gamma = 1/hum_t_inf
        self.mutation_p = mutation_p
        self.die_p = 1/mosq_t_inf
        self.K = K
        self.r = r
        self.population_hum = []
        self.population_mos_S = self.n_mosquitoes-self.init_inf_mos
        self.population_mos_I = self.init_inf_mos
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        # TODO convertir esto en un diccionario donde las llaves sean los genotipos y los valores el número de veces
        # que estan
        self.genotypes = pd.DataFrame(columns=["A", "C", "T", "G"])
        self.initialize_pop_hum()

    def initialize_pop_hum(self):
        # Based on the number of humans fills the list population_hum with objects of class human. Each one with
        # its id, its initial state (S) and an empty string as genotype.
        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state, "")
            self.population_hum.append(human)
        
        # From the list created above, it selects randomly the number of initial infected humans indicated in init_inf_hum
        # and in sample_infected are saved the humans picked randomly
        sample_infected = random.sample(self.population_hum, self.init_inf_hum)
        """ 
        For each human picked above, gets its index in the list population_hum, and then updates the state of this human
        in the list, from S to I
        """
        #TODO aqui ya no se infectan con un genotype the random choice sino del pool.
        
        for human in sample_infected:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            self.population_hum[index].genotype = random.choice(opt)
        """
        After updating the list population_hum with its corresponding states updates the dataframe data
        filling the ts as 0 for all humans and filling all other info (id, state, genotype) correspondingly.  
        """  
        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state, self.population_hum[i].genotype)
        
        # TODO REVISAR Creates conteos y genotype_counting with 0    
        """
        Runs the function conteos for the initial time step (0), as well as the function genotype_counting.
        """
        self.conteos(0)
        self.genotype_counting(0)
        
    def change_state(self):   
        """
        Updates the state of humans if necessary. 
        """

        
        
        """
        
        **Mosquito infection, human recovery and genotype mutation.**
        
        For each human in the list population_hum checks if is infected. 
        If it is, it is checked for each susceptible mosquito a random probability for encounter and if 
        this one is greater than encounter_p then it is checked a random probability of biting, if this one
        is greater than biting_p, the population of infected mosquitoes is updated with one infected mosquito more
        and the population of susceptible mosquitoes is updated with one susceptible mosquito less
        
        Then for each infected human checks if a random probability is greater than gamma, and if it is, the human
        is recovered and it is updated its state. If it is not, it checks if a random probability is greater than
        mutation_p, it it is, it updates the genotype of the human by the new mutated genotype.
        
        **Human infection.**
        
        For each human in the list population_hum checks if is susceptible. 
        If it is, it is checked for each infected mosquito a random probability for encounter and if 
        this one is greater than encounter_p then it is checked a random probability of biting, if this one
        is greater than biting_p, the state of the human is updated to infected and it is given a genotype to it.
        This genotype is the most frequent genotype of the population. Once that human it is infected, it passes
        to the next one.    
        
        **Mosquito dead**
        
        For each infected mosquito it checks if a random probability is greater than die_p, if it is, it 
        substracts 1 from the population_mos_I.
        
        """
        #TODO revisar con Mao esto porque implica que todas las personas infectadas están rodeadas de todos los
        # mosquitos susceptibles. Igual se tiene la probabilidad de encounter pero por si acaso.
        
        for human in self.population_hum:  
            if human.state == "I": 
                for mosquito in range(self.population_mos_S):    
                    if random.random() > self.encounter_p:                    
                        if random.random() > self.biting_p:
                            self.population_mos_I+=1
                            self.population_mos_S-=1
                               
                if random.random () > self.gamma:
                    human.state = "R"
                
                #TODO revisar esto, porque si sí muta, tiene que actualizarse el diccionario de los conteos de genotipos
                # o el data frame y quitar ese viejo. Además toca mirar como hacer que eso mute y que no mute un resto 
                #porque o sino al distancia genetica se jode.
                
                elif random.random () > self.mutation_p:
                    human.genotype = random.choice(opt)
                                    
            elif human.state == "S":
                for mosquito in range(self.population_mos_I):
                    if random.random() > self.encounter_p:
                        if random.random() > self.biting_p:
                            human.state ="I"
            
                            #TODO actualizar esto para que pues siga siendo el most frequent pero también mejor que lo saque de un diccionario o algo
                            # asi.
                            human.genotype = self.data["Genotype"].value_counts().idxmax()
                            # Si se infecta pues se asume que nadie más lo pica
                            break
                    
        for mosquito in range(self.population_mos_I):
            if random.random() > self.die_p:
                self.population_mos_I-=1
                    
    # Puede modificarlo teniendo transmisión vertical.
    
    def vital_dynamics(self):
        """
        Based on a logistic model for the mosquito population creates vital dynamics, this means that
        per time step find the new number of mosquitoes.
        """                
        # dN/dt = rN*(1-(N/k))
        N = self.population_mos_I + self.population_mos_S
        mosq = round((self.r*N)*(1-(N/self.K)))
        self.population_mos_S+=mosq  
                           
    def update(self, t):
        """
        In each step of time runs the function of changing the state of both mosquitoes and humans.
        Then, based on this, updates the dataframe data and calls the vital_dynamics function. 
        """
        self.change_state()
        for human in self.population_hum: 
            self.data.loc[len(self.data)] = (t,  human.id, human.state, human.genotype)   
        self.vital_dynamics()
    
    def conteos(self, t):
        """
        In each time step counts the number of humans in each of the states and updates
        the data frame counts. 
        """
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
                    
    #TODO falta mirar si dejar este conteo asi o en un diccionario o como
    
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
        """
        Corre todo el modelo en los diferentes pasos de tiempo. Con _init_ deja todas las condiciones iniciales, 
        correspondientes el tiempo 0. Así que run corre update, conteos y genotype_counting desde el ts 1. 
        Returns the data frame with the information of the agents by each time step, as well as the counts for 
        each state and for each genotype.
        """
        #mosquitos = [len(self.population_mos)]
        for t in range(1, n_steps):
            self.update(t)
            self.conteos(t)
            self.genotype_counting(t)
           # mosquitos.append(len(self.population_mos))
        return self.data, self.counts, self.genotypes

#Esto es pa ver cuanto se demora el código.

start_time = time.time()
# (self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, daysCured, mutation_p, K, r)
sims = 10
dias = 60
estados = 3
#x,y,z = simulaciones, tiempos, estado
matriz = np.zeros((sims, dias, estados))
for i in range(sims):
    model = SIRmodel(100, 100, 5, 10, 0.95, 0.95, 14, 0.2, 1500, 1/10, 4)
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


# Aca viene toda esa vaina de los intervalos de confianza

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

# Plotting
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






            
