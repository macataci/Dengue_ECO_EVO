from opqua.model import Model

my_model = Model()
my_model.newSetup('my_setup', preset='host-host')
my_model.newPopulation('my_population', 'my_setup', num_hosts=100)
my_model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )
my_model.run(0,100)
data = my_model.saveToDataFrame('../../output/my_model.csv')
graph = my_model.compartmentPlot('../../figs/my_model.png', data)

print(my_model.populations)
