import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
from Bio import Align
from statistics import mean 
import time
from collections import Counter
import statistics as stats

"""
En este script se presentarán las siguientes características:
- Se tendrá un modelo SIR en el cual los humanos se infectan al contacto con mosquitos infectados.
- Igualmente, los mosquitos se infectan al picar a humanos infectados.
- Los mosquitos serán cajas y no agentes dado a que no almacenaremos información relevante de ellos, tipo el genoma.
- El genoma asignado a los humanos saldrá de un pool creado inicialmente. Este pool corresponde a genotipos de un mismo serotipo, por lo que 
    sus distancias genéticas son pequeñas. El genoma asignado será el de mayor frecuencia. Ahora se manejan aminoacidos.
- Se tendrá mutación (pero pequeñita) y cuando esto suceda, esta nueva secuencia se añadirá al pool y se quitará la que se tenía antes de la mutación.
- También en caso de recuperación, ese genoma sale del pool, y se recuperan con una probabilidad uniforme. 
"""

#TODO revisar dinámicas vitales de humanos y mosquitos
#TODO mirar si crear una clase de tipo genotype 
#TODO definir escalas de tiempo

class Human:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
               
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, K, r, mosq_t_inf, amount, length, i_mut_region, f_mut_region, threshold):
        
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
                amount        (int): number of sequences in the initial pool
                length        (int): length of the genome
                i_mut_region  (int): starting position of the mutation region
                f_mut_region  (int): ending position of the mutation region
                threshold   (float): threshold that defines similarity
            """
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
        self.K = K
        self.r = r
        self.population_hum = []
        self.population_mos_S = self.n_mosquitoes-self.init_inf_mos
        self.population_mos_I = self.init_inf_mos
        
        # related to the genotype pool
        
        self.amount = amount
        self.length = length
        self.threshold = threshold
        self.i_mut_reg = i_mut_region
        self.f_mut_reg = f_mut_region
        self.genotype_counts = []
        self.dic_inititalpool = self.initial_pool()
        #self.letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        #related to humans
        
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        
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
        in the list, from S to I and assigns the secuence based on the frequencies of the initial pool, since these
        are the initial infected people
        """
        for human in sample_infected:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            self.population_hum[index].genotype = self.picking_from_pool(self.dic_inititalpool)
        """
        After updating the list population_hum with its corresponding states updates the dataframe data
        filling the ts as 0 for all humans and filling all other info (id, state, genotype) correspondingly.  
        """  
        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state, self.population_hum[i].genotype)
        """
        Runs the function conteos for the initial time step (0)
        """
        self.conteos_SIR(0)
        self.genotype_counting(0)
     
    # NO LA USO   
    def calculate_similarity (self, seq1, seq2):
        aligner = Align.PairwiseAligner(match_score=1.0)
        scoref = aligner.score(seq1, seq2)
        proportion = scoref/self.length
        return proportion
    
    def generate_random_string(self, length):
        letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        # Randomly choose characters from letters for the given length of the string
        random_string = ''.join(random.choice(letters) for i in range(length))
        return random_string 
      
    def initial_pool(self):
        # TODO solo se compara con la base, faltaría que se compararan entre todas pero ahi tiempo
        # computacional crece
        base = self.generate_random_string(self.length)
        pool = [base]
        for i in range (self.amount-1):
            sequence = self.generate_random_string(self.length)
            pool.append(sequence)
        conteo_t0 = Counter(pool)
        return conteo_t0
    
    def mutation(self, sequence, start, end):
        letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        input_seq = list(sequence)
        copy = input_seq.copy()
        posiciones = list(range(start,end+1))
        opciones = list(range(1,len(posiciones)+1))
        cuantos = random.choice(opciones)
        escogidos = random.sample(posiciones, cuantos)
        for i in escogidos:
            letter = random.choice(letters)
            while input_seq[i] == letter:
                letter = random.choice(letters)
            copy[i] = letter
        output_seq = ''.join([elem for elem in copy])
        return output_seq
        
    def picking_from_pool (self, dic):
        
        #TODO revisar esto bien, tipo qué hacer cuando ya todo el mundo es susceptible pero alguien se infecta
        if not dic: 
            #ant_ant = len(self.genotype_counts)-3
            dic =  self.initial_pool()
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
                    human.genotype = ""
                
                elif random.random () > self.mutation_p:
                    # TODO aqui va a comparar con cualquiera, falta que compare con todas, o que compare con la más
                    # repetida, y ahí uso la funcion de picking from the pool
                    new_seq = self.mutation(human.genotype, self.i_mut_reg, self.f_mut_reg)
                    #index_dic_anterior = len(self.genotype_counts)-1
                    #dic_anterior = self.genotype_counts[index_dic_anterior]
                    #seq_base = random.choice(list(dic_anterior.keys()))
                    #similarity = self.calculate_similarity(seq_base, new_seq)
                    #while similarity < self.threshold:
                     #   new_seq = self.mutation(human.genotype, self.i_mut_reg, self.f_mut_reg)
                      #  similarity = self.calculate_similarity(seq_base, new_seq)
                    human.genotype = new_seq
                                    
            elif human.state == "S":
                for mosquito in range(self.population_mos_I):
                    if random.random() > self.encounter_p:
                        if random.random() > self.biting_p:
                            human.state ="I"
                            index_dic_anterior = len(self.genotype_counts)-1
                            max_tant = self.picking_from_pool(self.genotype_counts[index_dic_anterior])
                            human.genotype = max_tant
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
                    
    #TODO falta mirar si dejar este conteo asi o en un diccionario o como
    
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
           # mosquitos.append(len(self.population_mos))
        return self.data, self.counts, self.genotype_counts


#Esto es pa ver cuanto se demora el código.

start_time = time.time()
# (n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, hum_t_inf, mutation_p, K, r, mosq_t_inf, amount, length, i_mut_reg, f_mut_reg, threshold):
#model = SIRmodel(100, 600, 5, 10, 0.9, 0.9, 2, 0.2, 1500, 1/10, 4, 10, 6, 3)
#df_data, df_conteos, list_dic_genotypes = model.run(10)
sims = 30
dias = 20
estados = 3
#x,y,z = simulaciones, tiempos, estado
probs = np.arange(0.1,1,0.1)
probs = np.append(probs, [0.95,0.97])
matriz = np.zeros((sims, dias, estados))

import os
dir_path = os.path.dirname(os.path.realpath(__file__))
ruta = os.path.join(dir_path, 'heatmap_vals_30sims.csv')
df_heatmaps = pd.DataFrame(columns=["Value_EB", "Value_R", "Genome_diversity"])
medianas = []
for j in probs:
    for k in probs:
        temporal = []
        for i in range(sims):
            union = []
            model = SIRmodel(500, 1500, 50, 150, j, j, k, 0.2, 3000, 1/10, 4, 20, 60, 2, 15, 0.6)
            df_data, df_conteos, list_dic_genotypes = model.run(dias)
            for l in range(len(list_dic_genotypes)):
                union += list(list_dic_genotypes[l].keys())
            union = np.unique(union)
            temporal.append(len(union))
        
        mediana = stats.median(temporal)
        df_heatmaps.loc[len(df_heatmaps)] = (j, k, mediana) 
            
df_heatmaps.to_csv(ruta, index=False) 
#print(df_heatmaps)

