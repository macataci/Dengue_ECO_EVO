from opqua.model import Model

my_model = Model()
my_model.newSetup('first_setup', preset='vector-borne')
my_model.newPopulation('humans1', 'first_setup', num_hosts=100, num_vectors=20)
my_model.addPathogensToHosts( 'my_population',{'agttgttagt':20}, {'ctacgtggac':30}, {'ctacgtggac':60} )
my_model.run(0,100)
data = my_model.saveToDataFrame('../../output/DENGUIX/my_model.csv')
graph = my_model.compartmentPlot('../../figs//DENGUIX/my_model.png', data)

print(my_model.populations)
