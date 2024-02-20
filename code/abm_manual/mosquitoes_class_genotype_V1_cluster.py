import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#from IPython.display import display
from Bio import Align
from statistics import mean 
import time
from collections import Counter
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


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
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, coinfection):
        
        self.n_humans = n_humans
        self.n_mosquitoes = n_mosquitoes
        self.init_inf_hum = init_inf_hum
        self.init_inf_mos = init_inf_mos
        self.encounter_p = encounter_p
        self.biting_p = biting_p
        #self.gamma = 1/hum_t_inf
        self.gamma = recovery_p
        self.mutation_p = mutation_p
        self.die_p = 1/mosq_t_inf
        #self.K = K
        #self.r = r
        self.population_hum = []
        self.population_mos = []
        

        
        self.amount = amount
        self.length = length
        self.genotype_counts = []
        self.coinfection = coinfection
        self.dic_initialpool = self.initial_pool()
        

        
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        
        self.initialize_pop_hum()
        self.initialize_pop_mos()

    def initialize_pop_hum(self):
        
        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state, "")
            self.population_hum.append(human)

        
        sample_infected_h = random.sample(self.population_hum, self.init_inf_hum)
        

        for human in sample_infected_h:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            self.population_hum[index].genotype = self.picking_from_pool(self.dic_initialpool)
            

        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state, self.population_hum[i].genotype)

        self.conteos_SIR(0)
        self.genotype_counting(0)
     
     
    def initialize_pop_mos(self):

        for i in range(self.n_mosquitoes):
            state = "S"
            mosquito = Mosquito(i, state, "")
            self.population_mos.append(mosquito)
        

        
        sample_infected_m = random.sample(self.population_mos, self.init_inf_mos)

        
        for mosquito in sample_infected_m:
            index = self.population_mos.index(mosquito)
            self.population_mos[index].state = "I"
            self.population_mos[index].genotype = self.picking_from_pool(self.dic_initialpool)
            
 
    
    def calculate_similarity (self, seq1, seq2):
        aligner = Align.PairwiseAligner(match_score=1.0)
        scoref = aligner.score(seq1, seq2)
        proportion = scoref/self.length
        return proportion
    
    def generate_random_string(self, length):
        letters = ["C", "G", "A", "U"]
        random_string = ''.join(random.choice(letters) for i in range(length))
        return random_string 
      
    def initial_pool(self):
        
        threshold = 0.85
        
        base = self.generate_random_string(self.length)
        pool = [base]
        for i in range (self.amount-1):
            sequence = self.generate_random_string(self.length)
            similarity = self.calculate_similarity(base, sequence)
            while similarity < threshold:
                sequence = self.generate_random_string(self.length)
                similarity = self.calculate_similarity(base, sequence)
            pool.append(sequence)
        conteo_t0 = Counter(pool)
        return conteo_t0
    
    def mutation(self, sequence, genome_length):
        letters = ["C", "G", "A", "U"]
        input_seq = list(sequence)
        copy = input_seq.copy()
        n, p = 1, 0.1 # number of trials, probability of each trial
        s = np.random.binomial(n, p, genome_length)
        for i in range(genome_length):
            if s[i]==1:
                letter = random.choice(letters)
                while input_seq[i] == letter:
                    letter = random.choice(letters)
                copy[i] = letter
        output_seq = ''.join([elem for elem in copy])
        return output_seq
                
    def picking_from_pool (self, dic):
        
        
        a = 2
        while not dic:
            a+=1
            dic = self.genotype_counts[len(self.genotype_count)-a]
        keys = dic.keys()
        values = dic.values()
        if max(values)==1:
            sequence = random.choice(list(keys))
        else:
            sequence_int = list({i for i in dic if dic[i]==max(values)})
            if len(sequence_int) == 1:
                sequence = sequence_int[0]
            else:
                sequence = random.choice(sequence_int)            
        return sequence
    
    def change_state(self):   
        
        for human in self.population_hum: 
            if human.state == "I": 
                for mosquito in self.population_mos:
                    if mosquito.state == "S":     
                        if random.random() > self.encounter_p:                    
                            if random.random() > self.biting_p:
                                mosquito.state = "I"
                                mosquito.genotype = human.genotype
                               
                if random.random () > self.gamma:
                    human.state = "R"
                    human.genotype = ""
                
                elif random.random () > self.mutation_p:                    
                    new_seq = self.mutation(human.genotype, self.length)
                    human.genotype = new_seq
                                    
            elif human.state == "S":
                for mosquito in self.population_mos:
                    if mosquito.state == "I":
                        if random.random() > self.encounter_p:
                            if random.random() > self.biting_p:
                                    human.state = "I"
                                    human.genotype = mosquito.genotype
                                    break

                                
        for mosquito in self.population_mos:
            if mosquito.state == "I":
                if random.random() > self.die_p:
                    self.population_mos.remove(mosquito)
                    


    """
    def vital_dynamics(self):
        
        Based on a logistic model for the mosquito population creates vital dynamics, this means that
        per time step find the new number of mosquitoes.
                     
        # dN/dt = rN*(1-(N/k))
        N = len(self.population_mos)
        mosq = round((self.r*N)*(1-(N/self.K)))
        self.population_mos_=mosq  
    """
                          
    def update(self, t):
        """
        In each step of time runs the function of changing the state of both mosquitoes and humans.
        Then, based on this, updates the dataframe data and calls the vital_dynamics function. 
        """
        self.change_state()
        for human in self.population_hum: 
            self.data.loc[len(self.data)] = (t,  human.id, human.state, human.genotype)   
        #self.vital_dynamics()
    
    def conteos_SIR(self, t):
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
                    
    
    def genotype_counting (self, t):
        
        genotypes_list = list(self.data.loc[(self.data["Time_step"]==t) & (self.data["State"]=="I")]["Genotype"])
        conteo_dict = Counter(genotypes_list)
        self.genotype_counts.append(conteo_dict)
        
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
            self.conteos_SIR(t)
            self.genotype_counting(t)
        return self.data, self.counts, self.genotype_counts



start_time = time.time() 
sims = 3
dias = 30
estados = 3
# Parametros modelo
n_humans = 500
n_mosquitoes = 1500
init_inf_hum = 50
init_inf_mos = 150
encounter_p = 0.97
biting_p = 0.97
recovery_p = 0.8
mutation_p = 0.2
#r = 1500
mosq_t_inf = 1/10
amount = 20
length = 20
coinfection = 0

matriz = np.zeros((sims, dias, estados))
for i in range(sims):        
    model = SIRmodel(n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, coinfection)
    df_data, df_conteos, list_dic_genotypes = model.run(dias)
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

# Plotting
colors = ["#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

figure, axis = plt.subplots(1,2)
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

union = []
for i in range(len(list_dic_genotypes)):
    union += list(list_dic_genotypes[i].keys())
union = np.unique(union)

df_genotypes = pd.DataFrame(columns=["Genotype", "Time", "Count", "Freq"])
for i in range(len(union)):
    for j in range(len(list_dic_genotypes)):
        if list_dic_genotypes[j].get(union[i]) is None:
            numero = 0
        else: 
            numero = list_dic_genotypes[j].get(union[i])
        total = sum(list_dic_genotypes[j].values())
        if total == 0:
            freq  = 0
        else: 
            freq = numero/total
        df_genotypes.loc[len(df_genotypes)] = (union[i], j, numero, freq) 
    freqs_seq = list(df_genotypes.loc[df_genotypes["Genotype"]==union[i]]["Freq"])
    promedio = np.mean(freqs_seq)
    if promedio >= 0.05:
        axis[1].plot(times, freqs_seq, label=union[i])
    else:
        axis[1].plot(times, freqs_seq)        
axis[1].set_xlabel('Time Step')
axis[1].set_ylabel('Genotype frequence')
axis[1].set_title('Genotype dynamics')
axis[1].legend(loc='upper right')
timef = time.time() - start_time
#print("--- %s seconds ---" % (time.time() - start_time))
with open("time.txt", "w") as file:
    # Write the text data to the file
    file.write(str(timef))
#plt.show()
plt.savefig('sample_plot.png')



