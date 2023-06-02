#mesa library

#Agentes - humans, mosquitos
#eventos: birth, death
#como entran los atributos?
#data collector para almacenar
#evento de me infecté o no? me recuperé o no?

#cantidad de gente en cada cajita
#bernoullli - binomial
#revisar distrib - #binomial de # de trials

#deltat = tasa*delta t me da alho similar a p
# poisson - # ebentos en un tiempo
# opqua tenia gillesies - metodo integracio

# mutacion
# revisar cosas biologicas
# wright-fischer, drift, mutacion y el otro
# poner mutacipn de genotipos en mis infectados, mirar el de jaime
#revisar oisson y revisar todo distribuciones
# revisar markov
# tengo un genoma y ese muta
###implementar markov


#100-500 agentes, algunos infectados
# y cada infectaado tiene que tener su genotipo de infeccion asociado
# la secuencia que tengo es proporcional a la cantidad que ella tiene en la poblacion.
# o sea tipo me infecto y lo que se me pega es al azar pero un poco ligado a proporción del total en la población



#evolu
#agentes - secuencias


#todo  estados: S, E, I, R
#ODE (aggregated) and ABM (disaggregated)

# Config
tiempo = 3600 # 3600 ticks = 3600 seconds = 1 hour
#I guess que parametros de infeccion, recuperacion, etc

#CLASS
from mesa import Agent

#todo revisar paquetes
class Human(Agent):
    def __init__(self, unique_id, model):
        #parametros ahora si r, i ?
        super().__init__(unique_id, model)
        # Time required to process the customer's transaction
        self.service_time = ss.poisson(45).rvs()
        # Time of arrival at queue
        self.entry_time = np.int(ss.beta(3, 3).rvs() * ticks) + 1
        self.balk_tolerance = ss.poisson(5).rvs() + 1
        # Whether or not the customer has arrived at the queue
        self._arrived = False
        self._chosen_counter = None
        self._q_entry = None
        self._q_exit = None
        # Start time when customer is being served
        self._service_entry = None
        # End time
        self._service_exit = None

    #metodos tipo: se encontró con infectados?
    # no se si las probs de infeccion y eso van como parmetros tipo en init o como metodos.
    # diria que parametros, pero entonces los metodos ya serían las cosas que tenía pablo tipo, mudó de poblacion
    # metodos - eventos

    def select_counter(self):
        self._arrived = True
        # Queue at shortest counter
        self._chosen_counter_idx = np.argmin([
            len(counter.queue) for counter in self.model.counters])
        self._chosen_counter = self.model.counters[
            self._chosen_counter_idx]
        # Balk if there are too many people at the counter
        if len(self._chosen_counter.queue) < self.balk_tolerance:
            self._chosen_counter.queue.append(self)
            self._q_entry = self.model._current_tick

    def pay_n_leave(self):
        self._service_exit = self.model._current_tick
        self._chosen_counter.active_customer = None

    #bueno esta aquí no sé

    def step(self):
        if (self._arrived == False) & \
                (self.model._current_tick >= self.entry_time):
            self.select_counter()
        elif isinstance(self._service_entry, int):
            if self.model._current_tick - self._service_entry \
                    == self.service_time:
                self.pay_n_leave()