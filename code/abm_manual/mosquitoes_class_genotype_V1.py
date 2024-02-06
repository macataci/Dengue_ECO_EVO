import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display
from Bio import Align
from statistics import mean 
import time
from collections import Counter


"""
import statistics
values = [6.28, 7.5, 7.77, 5.11, 7.8, 6.57]
mean = statistics.mean(values)
Median = statistics.median(values)
varianza = statistics.pvariance(values)
desv = statistics.pstdev(values)
print(mean, Median, varianza, desv)

"""

#TODO falta losgistic model
#evolution rates (reviar mosquito y human) - revisar espacio temporal 
#samplear de la mediana que tengo de esas evolutionary rates y hacer una distribucion
# mosquitos tienen eventos de nacimiento y muerte
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
    sus distancias genéticas son pequeñas. El genoma asignado será el de mayor frecuencia. Ahora se manejan aminoacidos.
    #TODO revisar si poner bases, aas o que
- Se tendrá mutación (pero pequeñita) y cuando esto suceda, esta nueva secuencia se añadirá al pool de ese tiempo y se quitará la que se tenía antes de la mutación.
- También en caso de recuperación, ese genoma sale del pool, y se recuperan con una probabilidad uniforme. 
"""

#TODO revisar dinámicas vitales de humanos y mosquitos
#TODO mirar si crear una clase de tipo genotype 
#TODO definir escalas de tiempo

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
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, i_mut_region, f_mut_region, coinfection):
        
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
                K              (int): carrying capacity mosquitoes  YA NO
                r            (float): constant of proportionality YA NO
                mosq_t_inf   (float): days that a mosquito is infected
                amount         (int): number of sequences in the initial pool
                length         (int): length of the genome
                i_mut_region   (int): starting position of the mutation region
                f_mut_region   (int): ending position of the mutation region
                coinfection   (bool): whether there is coinfection (1) or not(0)
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
        
        # related to the genotype pool
        
        self.amount = amount
        self.length = length
        self.i_mut_reg = i_mut_region
        self.f_mut_reg = f_mut_region
        self.genotype_counts = []
        self.coinfection = coinfection
        self.dic_initialpool = self.initial_pool()
        
        # related to humans
        
        self.data = pd.DataFrame(columns=["Time_step", "Id", "State", "Genotype"])
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
            human = Human(i, state, "")
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
            self.population_hum[index].genotype = self.picking_from_pool(self.dic_initialpool)
            
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
     
     
    def initialize_pop_mos(self):
        """ 
        Based on the number of mosquitoes fills the list population_mos with objects of class human. Each one with
        its id, its initial state (S) and an empty string as genotype.
        """ 
        
        for i in range(self.n_mosquitoes):
            state = "S"
            mosquito = Mosquito(i, state, "")
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
            self.population_mos[index].genotype = self.picking_from_pool(self.dic_initialpool)
            
        """
        Since we are not really interested on checking the genotype evolution or epidemic dynamics in mosquitoes, we do not
        save any kind of information on dataframes as we do with humans
        """    
    
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
        
        # TODO solo se compara con la secuencia base, faltaría que se compararan entre todas pero ahi tiempo
        # computacional crece
        
        # La diferencia maxima entre aa es de 3%, y para bases no excede el 6%
        # Por lo que deben tener una similaridad mínima de 97% o 94% respectivamente
        threshold = 0.6
        
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
    
    def mutation(self, sequence, start, end):
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
                for mosquito in range(self.population_mos):
                    if mosquito.state == "S":     
                        if random.random() > self.encounter_p:                    
                            if random.random() > self.biting_p:
                                mosquito.state = "I"
                                mosquito.genotype = human.genotype
                               
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
                for mosquito in range(self.population_mos):
                    if mosquito == "I":
                        if random.random() > self.encounter_p:
                            if random.random() > self.biting_p:
                                human.state ="I"
                                human.genotype = mosquito.genotype
                                """TODO 
                                if self.coinfection == 0:
                                    break
                                else:
                                """
                                
        for mosquito in range(self.population_mos):
            if mosquito == "I":
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
        return self.data, self.counts, self.genotype_counts


#Esto es pa ver cuanto se demora el código.

start_time = time.time() 
sims = 1
dias = 20
estados = 3
# Parametros modelo
n_humans = 100
n_mosquitoes = 200
init_inf_hum = 1
init_inf_mos = 2
encounter_p = 0.9
biting_p = 0.9
recovery_p = 0.2
mutation_p = 0.2
#r = 1500
mosq_t_inf = 1/10
amount = 5
length = 20
i_mut_region = 4
f_mut_region = 10
coinfection = 0

#x,y,z = simulaciones, tiempos, estado
matriz = np.zeros((sims, dias, estados))
for i in range(sims):
    model = SIRmodel(n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, recovery_p, mutation_p, mosq_t_inf, amount, length, i_mut_region, f_mut_region, coinfection)
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
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()




"""
Notas:
Si cuando muta lo pongo a comparar, y la region es muy pequeña, y solo tengo 1 SNP por mutacion por ts, puede que se quede en un loop 
infinito porque nunca pase el umbral, incluso teniendo como region de mutacion toda la sequencia.
"""


