from opqua.model import Model

my_model = Model()
my_model.newSetup('first_setup', preset='vector-borne', protection_upon_recovery_host=[0,3], possible_alleles='acgt',
                  mutate_in_host=2.06e-7)
my_model.newPopulation('pop1', 'first_setup', num_hosts=100, num_vectors=1000)
#subir mosquitos
#10% personas vs mosquitos
my_model.addPathogensToVectors('pop1',{'agttgttagt':100, 'ttctaacagt':50, 'tcgaccgtct':50, 'ggcgaagaga':50})
my_model.newIntervention(240, 'addPathogensToVectors', [ 'pop1', {'ttttatagca':5, 'agttgttagt':20}])
my_model.newIntervention(400, 'addPathogensToVectors', [ 'pop1', {'gggctcattc':5}])

my_model.newSetup('second_setup', preset='vector-borne', protection_upon_recovery_host=[0,3], possible_alleles='acgt',
                  mutate_in_host=2.06e-1)

my_model.newIntervention( 423, 'setSetup', ['pop1', 'second_setup'] )
my_model.run(0,730)
data = my_model.saveToDataFrame('../../../output/DENGUIX/my_model.csv')
#patho_data = my_model.getCompositionData(data=data, save_data_to_file='../../../output/DENGUIX/pathogens.csv', vectors=True, track_specific_sequences=['agttgttagt'])
pop_plot_I =my_model.populationsPlot( '../../../figs//DENGUIX/pop_plot_I.png',data, vectors=True, y_label="Infected hosts and vectors", compartment="Infected")
pop_plot_R =my_model.populationsPlot( '../../../figs//DENGUIX/pop_plot_R.png',data, vectors=True, y_label="Recovered hosts and vectors", compartment="Recovered")
graph = my_model.compartmentPlot('../../../figs//DENGUIX/my_model.png', data, vectors=True,y_label="Hosts and vectors" )
composition_path = my_model.compositionPlot( '../../../figs/DENGUIX/composition_path.png', data,
                                             track_specific_sequences=['agttgttagt'], type_of_composition='Pathogens',
                                             vectors=True,y_label="Hosts and vectors" )
composition_prot=my_model.compositionPlot( '../../../figs/DENGUIX/composition_prot.png', data,
                                           track_specific_sequences=['agttgttagt'], type_of_composition='Protection',
                                           vectors=True,y_label="Hosts and vectors" )


graph_clustermap = my_model.clustermap(
    '../../../figs/DENGUIX/clustermap.png', data,
    save_data_to_file='../../../output/DENGUIX/pairwise_distances.csv',
    )
