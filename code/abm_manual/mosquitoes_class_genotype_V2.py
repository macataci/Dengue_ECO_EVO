import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#from IPython.display import display
from Bio import Align
from statistics import mean 
import time
from collections import Counter
import math


#TODO falta losgistic model
#evolution rates (reviar mosquito y human) - revisar espacio temporal 
#samplear de la mediana que tengo de esas evolutionary rates y hacer una distribucion
#mosquitos tienen eventos de nacimiento y muerte
#grafica distancia gentica y antigenica
#meterle 2 años a las simulaciones
#revisar inmunidades y ADE
#revisar coinfeccion
"""
En este script se presentarán las siguientes características:
- Se tendrá un modelo SIR en el cual los humanos se infectan al contacto con mosquitos infectados.
- Igualmente, los mosquitos se infectan al picar a humanos infectados.
- Los mosquitos serán tratados como agentes y se almacenará información como el id, estado de infección y genotipo
- El genoma asignado a los humanos y mosquitos saldrá de un pool creado inicialmente. Este pool corresponde a genotipos de un mismo serotipo, por lo que 
    sus distancias genéticas son pequeñas. El genoma asignado será el de mayor frecuencia. Ahora se manejarán nucleotidos y se hará la traduccion.
- Se tendrá mutación (pero pequeñita) y cuando esto suceda, esta nueva secuencia se añadirá al pool de ese tiempo y se quitará la que se tenía antes de la mutación.
- También en caso de recuperación, ese genoma sale del pool, y se recuperan con una probabilidad uniforme. 
- Se tiene en cuenta traduccion
- Se expande a mas de un serotipo
"""

#TODO revisar dinámicas vitales de humanos y mosquitos
#TODO mirar si crear una clase de tipo genotype 
#TODO definir escalas de tiempo
#TODO mi mosquito necesita almacenar su proteina?

class Mosquito:
    def __init__(self, id, state, serotype, genotype): 
        self.id = id
        self.state = state
        self.serotype = serotype
        self.genotype = genotype
        
class Human:
    def __init__(self, id, state, serotype, genotype, protein): 
        self.id = id
        self.state = state
        self.serotype = serotype
        self.genotype = genotype
        self.protein = protein
             
class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, coinfection, serotypes):
        
        """
            Args:
                n_humans       (int): human population
                n_mosquitoes   (int): mosquito population
                init_inf_hum   (int): initial number of infected humans          
                init_inf_mos   (int): initial number of infected mosquitoes
                encounter_p  (float): probability of encounter of a human with a mosquito and viceverse
                bitting_p    (float): probability that the encounter results in a bite from the mosquito
                recovery_p    float): probablity of recovery
                mutation_p   (float): mutation rate of the genome
                mosq_t_inf   (float): days that a mosquito is infected
                amount         (int): number of sequences in the initial pool
                length         (int): length of the genome
                coinfection   (bool): whether there is coinfection (1) or not(0)
                seroypes       (int): number of serotypes
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
        #self.K = K
        #self.r = r
        self.population_hum = []
        self.population_mos = []
        
        # related to the genotypes/serotypes
        
        self.amount = amount
        self.length = length
        self.serotypes = serotypes
        self.genotype_counts = []
        self.protein_counts = []
        self.serotype_counts = []
        self.coinfection = coinfection
        #self.dic_initialpool = self.initial_pool()
        self.pools_list = self.serotype_creation()
        # related to humans
        
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Serotype", "Genotype", "Protein"])
        self.counts = pd.DataFrame(columns=["S", "I", "R"])
        
        self.initialize_pop_hum()
        self.initialize_pop_mos()

    def initialize_pop_hum(self):
        """
        Based on the number of humans fills the list population_hum with objects of class human. Each one with
        its id, its initial state (S) and an empty string as genotype.
        """
        
        for i in range(self.n_humans):
            state = "S"
            human = Human(i, state, "", "", "")
            self.population_hum.append(human)
        
        """
        From the list created above, it selects randomly the initial infected humans, based on init_inf_hum
        and in sample_infected are saved the humans picked randomly
        """ 
        
        sample_infected_h = random.sample(self.population_hum, self.init_inf_hum)
        
        
        """
        For each human picked above, gets its index in the list population_hum, and then updates the state of this human
        in the list, from S to I and assigns the secuence based on the frequencies of the initial pool, since these
        are the initial infected people
        """
        
        for human in sample_infected_h:
            index = self.population_hum.index(human)
            self.population_hum[index].state = "I"
            selected_serotype = random.randint(0, self.serotypes-1)
            self.population_hum[index].serotype = selected_serotype + 1
            self.population_hum[index].genotype = self.picking_from_pool(self.pools_list[selected_serotype])
            self.population_hum[index].protein = self.translating_RNA_prot(human.genotype)
        
        """
        After updating the list population_hum with its corresponding states updates the dataframe data
        filling the ts as 0 for all humans and filling all other info (id, state, genotype) correspondingly.  
        """  
        
        for i in range(len(self.population_hum)):
            self.data.loc[i]=(0, self.population_hum[i].id, self.population_hum[i].state, self.population_hum[i].serotype, self.population_hum[i].genotype, self.population_hum[i].protein)
        """
        Runs the function conteos for the initial time step (0)
        """
        
        self.conteos_SIR(0)
        self.genotype_counting(0)
        self.protein_counting(0)
        self.serotype_counting(0)
     
    def initialize_pop_mos(self):
        """ 
        Based on the number of mosquitoes fills the list population_mos with objects of class human. Each one with
        its id, its initial state (S) and an empty string as genotype.
        """
        for i in range(self.n_mosquitoes):
            state = "S"
            mosquito = Mosquito(i, state, "", "")
            self.population_mos.append(mosquito)
        
        """ 
        From the list created above, it selects randomly the initial infected mosquitoes, based on init_inf_mos
        and in sample_infected are saved the mosquitoes picked randomly
        """ 
        
        sample_infected_m = random.sample(self.population_mos, self.init_inf_mos)
        
        """ 
        For each human picked above, gets its index in the list population_mos, and then updates the state of this mosquito
        in the list, from S to I and assigns the secuence based on the frequencies of the initial pool, since these
        are the initial infected mosquitoes
        """
        
        for mosquito in sample_infected_m:
            index = self.population_mos.index(mosquito)
            self.population_mos[index].state = "I"
            selected_serotype = random.randint(0, self.serotypes-1)
            self.population_mos[index].serotype = selected_serotype + 1
            self.population_mos[index].genotype = self.picking_from_pool(self.pools_list[selected_serotype])
            
        """
        Since we are not really interested on checking the genotype evolution or epidemic dynamics in mosquitoes, we do not
        save any kind of information on dataframes as we do with humans
        """    
    
    def translating_RNA_prot(self, rna_seq):
        # RNA codon table
        rna_codon = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
           "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
           "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
           "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
           "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
           "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }

        protein_string = ""
        # Generate protein string
        for i in range(0, len(rna_seq)-(3+len(rna_seq)%3), 3):
            if rna_codon[rna_seq[i:i+3]] == "STOP" :
                break
            protein_string += rna_codon[rna_seq[i:i+3]]

        # Print the protein string
        return protein_string
    
    def calculate_similarity (self, seq1, seq2):
        aligner = Align.PairwiseAligner(match_score=1.0)
        scoref = aligner.score(seq1, seq2)
        proportion = scoref/self.length
        return proportion
    
    def generate_random_string(self, length):
        #letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        letters = ["C", "G", "A", "U"]
        # Randomly choose characters from letters for the given length of the string
        random_string = ''.join(random.choice(letters) for i in range(length))
        return random_string 
    
    def find_differences(self, str1, str2):
        differences = []
        for i, (char1, char2) in enumerate(zip(str1, str2)):
            if char1 != char2:
                differences.append(i)
        return differences
    
    def find_similarities(self, str1, str2):
        similarities = []
        for i, (char1, char2) in enumerate(zip(str1, str2)):
            if char1 == char2:
                similarities.append(i)
        return similarities
    
    def serotype_creation (self):
        serotype_number = self.serotypes
        pools = []
        letters = ["C", "G", "A", "U"]
        base = self.generate_random_string(self.length)
        base_seqs = [base]
        boundaries = [0.65, 0.75]
        max_needed = (boundaries[1]*self.length)
        min_needed = (boundaries[0]*self.length)
        if serotype_number >=2:
            for i in range(serotype_number):
                sequence = self.generate_random_string(self.length)
                similarity = self.calculate_similarity(base, sequence)
                actual = (similarity*self.length)
                while similarity > boundaries[1]:
                    # Necesito reducir similaridad
                    needing_min = math.ceil(actual-min_needed)
                    needing_max = math.ceil(actual-max_needed)
                    to_change = random.randint(needing_min, needing_max)
                    which = random.sample(self.find_similarities(base, sequence), to_change)
                    sequence_chars = list(sequence)
                    for pos in which:
                        letra = random.choice(letters)
                        while letra != base[pos]:
                            sequence_chars[pos] = letra
                    sequence = ''.join(sequence_chars)
                    similarity = self.calculate_similarity(base, sequence)
                while similarity < boundaries[0]:
                    # Necesito aumentar similaridad
                    needing_min = math.ceil(min_needed- actual)
                    needing_max = math.ceil(max_needed - actual)
                    to_change = random.randint(needing_min, needing_max)
                    which = random.sample(self.find_differences(base, sequence), to_change)
                    sequence_chars = list(sequence)
                    for pos in which:
                        sequence_chars[pos] = base[pos]
                    sequence = ''.join(sequence_chars)
                    similarity = self.calculate_similarity(base, sequence)                   
                base_seqs.append(sequence)
        
        for i in range(serotype_number):
            pool = self.initial_pool(base_seqs[i])
            pools.append(pool)           
        return (pools)
    
    def initial_pool(self, base):
        # TODO solo se compara con la secuencia base, faltaría que se compararan entre todas pero ahi tiempo
        # computacional crece
        # La diferencia maxima entre aa es de 3%, y para bases no excede el 6%
        # Por lo que deben tener una similaridad mínima de 97% o 94% respectivamente
        # TODO quizás sea bueno poner codon de parada al final de todas
        threshold = 0.94
        pool = [base]
        for i in range (amount-1):
            sequence = self.generate_random_string(self.length)
            similarity = self.calculate_similarity(base, sequence)
            while similarity < threshold:
                min_needed = (threshold*self.length)
                actual = (similarity*self.length)
                needing_min = math.ceil(min_needed-actual)
                needing_max = math.ceil(self.length-actual)-1
                to_change = random.randint(needing_min, needing_max)
                which = random.sample(self.find_differences(base, sequence), to_change)
                sequence_chars = list(sequence)
                for pos in which:
                    sequence_chars[pos] = base[pos]
                sequence = ''.join(sequence_chars)
                similarity = self.calculate_similarity(base, sequence)
            pool.append(sequence)
        conteo_t0 = Counter(pool)
        return conteo_t0
    
    def mutation(self, sequence, genome_length):
        """
        Binomial para el largo de la secuencia y con eso miro si muta o no cada site
        
        
        Poisson y miro para cada valor entre 1 y el largo del genoma, y escojo el de mayor p
        Bootstrap y para cada site saco un valor aleaotrio de la distribucion. con lambda igual a la mediana de los rates de alguna distribucion
        Normal donde asigno los valores para cada 
        
        
        Median evolutionary rates by seroytpe (X10-4 subs/site/year):
        DENV1 = 6.57
        DENV2 = 6.98
        DENV3 = 8.55
        DENV4 = 7.91
        """
        #letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
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
        
        #TODO revisar esto bien, tipo qué hacer cuando ya todo el mundo es susceptible pero alguien se infecta.
        #SOLVED: se va devolviendo de a uno
        
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
        """
        Updates the state of humans if necessary. 
        """
        
        """
        **Mosquito infection, human recovery and genotype mutation**
        
        TODO revisar si mejor primer loop de humanos o s de mosquitos
        
        For each human in the list population_hum checks if is infected. 
        If it is, it is checked for each susceptible mosquito a random probability for encounter and if 
        this one is greater than encounter_p then it is checked a random probability of biting, if this one
        is greater than biting_p, the state of infection of the mosquito is updated to "I" and it is updated also
        the genotype value corresponding to the genotype of infection associated to the human that was bitten.
        
        Then for each infected human checks if a random probability is greater than gamma, and if it is, the human
        is recovered and it is updated its state. If it is not, TODO use evolutionary rates to introduce mutation
        
        **Human infection**
        
        For each human in the list population_hum checks if is susceptible. 
        If it is, it is checked for each infected mosquito a random probability for encounter and if 
        this one is greater than encounter_p then it is checked a random probability of biting, if this one
        is greater than biting_p, the state of the human is updated to infected and its associated genotype is the one
        associated to the infected mosquito. Since the human was already infected, it passes to the next human.
        TODO cambiar si coinfeccion
        
        **Mosquito dead**
        
        TODO mirar eventos de nacimiento y muerte
        For each infected mosquito it checks if a random probability is greater than die_p, if it is, it 
        delets this mosquito from the pupulation_mos.
        
        """
        #TODO revisar con Mao esto porque implica que todas las personas infectadas están rodeadas de todos los
        # mosquitos susceptibles. Igual se tiene la probabilidad de encounter pero por si acaso.
        
        for human in self.population_hum: 
            #cambiar esto a donde sea indexado 
            if human.state == "I": 
                for mosquito in self.population_mos:
                    if mosquito.state == "S":     
                        if random.random() > self.encounter_p:                    
                            if random.random() > self.biting_p:
                                mosquito.state = "I"
                                mosquito.serotype = human.serotype
                                mosquito.genotype = human.genotype
                               
                if random.random () > self.gamma:
                    human.state = "R"
                    human.serotype = ""
                    human.genotype = ""
                    human.protein = ""
                
                elif random.random () > self.mutation_p:
                    # TODO aqui va a comparar con cualquiera, falta que compare con todas, o que compare con la más
                    # repetida, y ahí uso la funcion de picking from the pool
                    # TODO asumiré que sigue siendo el mismo serotipo
                    
                    new_seq = self.mutation(human.genotype, self.length)
                    human.genotype = new_seq
                    human.protein = self.translating_RNA_prot(human.genotype)
                                    
            elif human.state == "S":
                for mosquito in self.population_mos:
                    if mosquito.state == "I":
                        if random.random() > self.encounter_p:
                            if random.random() > self.biting_p:
                                    human.state = "I"
                                    human.serotype = mosquito.serotype
                                    human.genotype = mosquito.genotype
                                    human.protein = self.translating_RNA_prot(human.genotype)
                                    break
                                    """
                                    TODO 
                                    if self.coinfection == 0:
                                         break
                                    else:
                                    """
                                
        for mosquito in self.population_mos:
            if mosquito.state == "I":
                if random.random() > self.die_p:
                    self.population_mos.remove(mosquito)
                    
    # Puede modificarlo teniendo transmisión vertical.
    
    
    #TODO voy a quitar logistico mientras, revisar
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
            self.data.loc[len(self.data)] = (t,  human.id, human.state, human.serotype, human.genotype, human.protein)   
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
                    
    #TODO falta mirar si dejar este conteo asi o en un diccionario o como
    
    def genotype_counting (self, t):
        
        genotypes_list = list(self.data.loc[(self.data["Time_step"]==t) & (self.data["State"]=="I")]["Genotype"])
        conteo_dict = Counter(genotypes_list)
        self.genotype_counts.append(conteo_dict)
        
    def protein_counting (self, t):
        
        proteins_list = list(self.data.loc[(self.data["Time_step"]==t) & (self.data["State"]=="I")]["Protein"])
        conteo_dict = Counter(proteins_list)
        self.protein_counts.append(conteo_dict)    
            
    def serotype_counting (self, t):
        
        serotype_list = list(self.data.loc[(self.data["Time_step"]==t) & (self.data["State"]=="I")]["Serotype"])
        conteo_dict = Counter(serotype_list)
        self.serotype_counts.append(conteo_dict)    
                   
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
            self.serotype_counting(t)
            self.genotype_counting(t)
            self.protein_counting(t)
        return self.data, self.counts, self.serotype_counts, self.genotype_counts, self.protein_counts


#Esto es pa ver cuanto se demora el código.

start_time = time.time() 
sims = 3
dias = 150
estados = 3
# Parametros modelo
n_humans = 50
n_mosquitoes = 150
init_inf_hum = 5
init_inf_mos = 15
encounter_p = 0.97
biting_p = 0.97
recovery_p = 0.8
mutation_p = 0.2
#r = 1500
mosq_t_inf = 1/10
amount = 20
length = 20
coinfection = 0
serotipos = 4

#x,y,z = simulaciones, tiempos, estado
matriz = np.zeros((sims, dias, estados))
for i in range(sims):        
    model = SIRmodel(n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, coinfection, serotipos)
    df_data, df_conteos, list_dic_serotypes, list_dic_genotypes, list_dic_proteins = model.run(dias)
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

figure, axis = plt.subplots(2,2, figsize=(20, 25))
times = range(dias)

# PLOT EPI
axis[0,0].plot(times, mediaS_total, label="S", color=colors[0])
axis[0,0].fill_between(times, statistics_95["LowS"], statistics_95["HighS"], alpha=0.2, color=colors[0])
axis[0,0].fill_between(times, statistics_50["LowS"], statistics_50["HighS"], alpha=0.3, color=colors[0])

axis[0,0].plot(times, mediaI_total, label="I", color=colors[1])
axis[0,0].fill_between(times, statistics_95["LowI"], statistics_95["HighI"], alpha=0.2, color=colors[1])
axis[0,0].fill_between(times, statistics_50["LowI"], statistics_50["HighI"], alpha=0.3, color=colors[1])

axis[0,0].plot(times, mediaR_total, label="R", color=colors[2])
axis[0,0].fill_between(times, statistics_95["LowR"], statistics_95["HighR"], alpha=0.2, color=colors[2])
axis[0,0].fill_between(times, statistics_50["LowR"], statistics_50["HighR"], alpha=0.3, color=colors[2])

#axis[0,0].set_xlabel('Time Step')
axis[0,0].set_ylabel('Number of Humans')
axis[0,0].set_title('Humans dynamics')
axis[0,0].legend()

# PLOT GENOTIPOS
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
        axis[0,1].plot(times, freqs_seq, label=union[i])
    else:
        axis[0,1].plot(times, freqs_seq)        
#axis[0, 1].set_xlabel('Time Step')
axis[0, 1].set_ylabel('Genotype frequence')
axis[0, 1].set_title('Genotype dynamics')
axis[0, 1].legend(loc='upper right')

# PLOT PROTEINAS

union_prot = []
for i in range(len(list_dic_proteins)):
    union_prot += list(list_dic_proteins[i].keys())
union_prot = np.unique(union_prot)

df_proteins = pd.DataFrame(columns=["Protein", "Time", "Count", "Freq"])
for i in range(len(union_prot)):
    for j in range(len(list_dic_proteins)):
        if list_dic_proteins[j].get(union_prot[i]) is None:
            numero = 0
        else: 
            numero = list_dic_proteins[j].get(union_prot[i])
        total = sum(list_dic_proteins[j].values())
        if total == 0:
            freq  = 0
        else: 
            freq = numero/total
        df_proteins.loc[len(df_proteins)] = (union_prot[i], j, numero, freq) 
    freqs_prot = list(df_proteins.loc[df_proteins["Protein"]==union_prot[i]]["Freq"])
    promedio = np.mean(freqs_prot)
    if promedio >= 0.05:
        axis[1,0].plot(times, freqs_prot, label=union_prot[i])
    else:
        axis[1,0].plot(times, freqs_prot)        
axis[1,0].set_xlabel('Time Step')
axis[1,0].set_ylabel('Protein frequence')
axis[1,0].set_title('Protein dynamics')
axis[1,0].legend(loc='upper right')


# PLOT SEROTIPOS
union_ser = []
for i in range(len(list_dic_serotypes)):
    union_ser += list(list_dic_serotypes[i].keys())
union_ser = np.unique(union_ser)
df_serotypes = pd.DataFrame(columns=["Serotype", "Time", "Count", "Freq"])

for i in range(len(union_ser)):
    for j in range(len(list_dic_serotypes)):
        if list_dic_serotypes[j].get(union_ser[i]) is None:
            numero = 0
        else: 
            numero = list_dic_serotypes[j].get(union_ser[i])
        total = sum(list_dic_serotypes[j].values())
        if total == 0:
            freq  = 0
        else: 
            freq = numero/total
        df_serotypes.loc[len(df_serotypes)] = (union_ser[i], j, numero, freq) 
    freqs_ser = list(df_serotypes.loc[df_serotypes["Serotype"]==union_ser[i]]["Freq"])
    promedio = np.mean(freqs_ser)
    if promedio >= 0:
        axis[1,1].plot(times, freqs_ser, label=union_ser[i])
    else:
        axis[1,1].plot(times, freqs_ser)        
axis[1, 1].set_xlabel('Time Step')
axis[1, 1].set_ylabel('Serotype frequence')
axis[1, 1].set_title('Serotype dynamics')
axis[1, 1].legend(loc='upper right')

timef = time.time() - start_time
print("--- %s seconds ---" % (time.time() - start_time))
#with open("time.txt", "w") as file:
    # Write the text data to the file
    #file.write(timef)
with open("prot.txt", "w") as file:
    # Write the text data to the file
    file.write(str(list_dic_proteins))    
with open("gen.txt", "w") as file:
    # Write the text data to the file
    file.write(str(list_dic_genotypes))    
    

plt.show()
#plt.savefig('sample_plot.png')