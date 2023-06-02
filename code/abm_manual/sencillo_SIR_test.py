#Enter Python code here and hit the Run button.
import random
import pandas as pd
from IPython.display import display

class Human:
    def __init__(self, id, state): 
        self.id = id
        self.state = state
# Init

n_humans = 5
init_inf_hum = 2
beta = 0.6
gamma = 0.3
population = []

data = pd.DataFrame(columns=["Time_step", "Id", "State"])
counts = pd.DataFrame(columns=["S", "I", "R"])

# Initialize
def initialize():
    for i in range(n_humans):
        state = "S"
        human = Human(i, state)
        population.append(human)
    
    sample_infected = random.sample(population, init_inf_hum)
    
    for human in sample_infected:
        index = population.index(human)
        population[index].state = "I"
    
    for i in range(len(population)):
        data.loc[i]=(0, population[i].id, population[i].state)
        
    conteos(0)
    print(data)
    print(counts)

#DESDE ACA PASA EN EL TIEMPO 1

# Change state
def change():
    for human in population:
        #Infection
        if human.state=="S":
            if random.random()< beta:
                human.state="I"
                
        #Recovery
        elif human.state=="I":
            if random.random()< gamma:
                human.state="R" 
def update(t):
    change()
    for human in population: 
        data.loc[len(data)] = (t,  human.id, human.state)

def conteos(t):
    conteo = data.where(data["Time_step"]==t)["State"].value_counts()
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
        
    counts.loc[t] = (S, I, R)     
    
def run(n_steps):
    initialize()
    for t in range(1, n_steps):
        update(t)
        conteos(t)
    display(data)
    display(counts)
    
run(10)
