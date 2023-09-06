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


class Human:
    def __init__(self, id, state, genotype): 
        self.id = id
        self.state = state
        self.genotype = genotype
        
        
        

class SIRmodel:
    def __init__(self, n_humans, n_mosquitoes, init_inf_hum, init_inf_mos, encounter_p,  biting_p, hum_t_inf, mutation_p, K, r, mosq_t_inf):
            """
            Args:
                n_humans       (int): human population
                n_mosquitoes   (int): mosquito population
                init_inf_hum   (int): initial number of infected humans          
                init_inf_mos   (int): initial number of infected mosquitoes
                encounter_p  (float): probability of encounter of a human with a mosquitoe and viceverse
                bitting_p    (float): probability that the encounter results in a bite
                
                

            """
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

