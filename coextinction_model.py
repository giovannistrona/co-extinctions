from scipy import trapz
from random import random,randrange,sample,shuffle
from numpy import arange,array,where,nan_to_num,zeros,mean,linspace
from copy import deepcopy
from scipy.spatial.distance import euclidean
from igraph import Graph
from numpy.random import normal,lognormal
from difflib import SequenceMatcher
from load import csv_

import string

###########################################################################################
### LOAD REAL WORLD DATASETS ON FOOD WEB STRUCTURE, CLIMATE, AND SPECIES THERMAL TOLERANCE
###########################################################################################

# 1) trophic level, specialization
tl_file=[map(float,i[:-1]) for i in csv_('calibration.csv')] 
# 2) climate data: lat,lon,mean T,min T,max T
env_file=[[float(i[0]),float(i[1]),float(i[3])/10.0,float(i[8])/10.0,float(i[7])/10.0] for i in csv_('./thermal_tolerance/bioclim_terrestrial.csv')[1:] if i[3]!='na'] 
# 3) plant thermal tolerance: maxT,minT
plant_file=[[float(i[3])/10.0,float(i[2])/10.0] for i in csv_('./thermal_tolerance/plants.csv')] 
# 4) endotherm thermal tolerance: tmaxT,minT
endo_file=[map(float,i[1:]) for i in csv_('./thermal_tolerance/endotherm.csv')] 
# 5) ectotherm thermal tolerance: maxT,minT
ecto_file=[map(float,i[1:]) for i in csv_('./thermal_tolerance/ectotherm_terrestrial.csv')] 

###########################################################################################



############################################
#IMPLEMENTATION OF SPECIES FUNCTIONAL TRAITS
############################################

'''
1) Generate an adjacency matrix for functional traits (each trait is one of 26 letters);
the value in cell ij indicates the the degree having the i-th trait enables a consumer 
to use a resource having the j-th trait. Positive values indicate that the trait i-th makes 
consumption of a resource having trait j easier for the consumer. Negative values indicate 
that the trait j constitutes a protection against trait i. 
Organisms will be assigned a string of 10 letters, ideally a proxy for their phenotype. 
The ability of a consumer to use a potential resource will be then quantified by summing up all 
the ij entries in the adjacency matrix for each i trait of the consumer and each j trait of the 
potential resource.
'''

ft = string.ascii_lowercase
ft_mat = zeros([26,26])
for i in range(26):
	for j in range(26):
		if i!=j:
			sign = sample([-1,1],1)[0]
			ft_mat[i][j] = sign*random()
		else:
			ft_mat[i][j] = 0


'''
2) Since we need to rescale the values stemming from the comparison of two phenotypes between 
0 and 1, 1M pairs of random phenotypes are generated and evaluated for resource-consumer 
compatibility. The minimum and maximum empirical compatibility (which depend on the trait 
adjacency matrix) are used in the rest of the analyses to rescale compatibility values.
'''
res = []
for rep in range(1000000):
	resource = sample(ft,randrange(1,10))			
	consumer = sample(ft,randrange(1,10))
	sc = 0
	for c_trait in consumer:
		for r_trait in resource:
			sc += ft_mat[ft.index(c_trait)][ft.index(r_trait)]#
	res.append(sc)


min_tc, max_tc  = min(res),max(res)



'''	
3) Define the function to evaluate a consumer's ability to use a resource, based on the
resource and consumer phenotypes (i.e. strings including a random sample of letters of size varying
between 1 and 10), the trait adjacency matrix, and the minimum/maximum compatibility values 
estimated empirically in the previous step. As explained above, compatibility is obtained by 
summing up all the ij entries in the  adjacency matrix for each i trait of the consumer and each j 
trait of the potential resource.
'''

def compute_tc(resource,consumer,ft,ft_mat, min_tc, max_tc):
	sc = 0
	for c_trait in consumer:
		for r_trait in resource:
			sc += ft_mat[ft.index(c_trait)][ft.index(r_trait)]
	val = 1-(max_tc-sc)/(max_tc-min_tc)
	if val>1:
		val = 1
	elif val<0:
		val = 0	
	return val

##########################################################################################
###DEFINE STRING MATCHING FUNCTION TO COMPARE PHENOTYPES
##########################################################################################

def phen_similarity(phen1,phen2):
	s = SequenceMatcher()
	s.set_seqs(phen1, phen2)
	return s.ratio()
	



##########################################################################################
###GENERATE GLOBAL SPECIES POOL
##########################################################################################
'''Each entry in the species pool is a list representing an organism's feature, including: 
unique identifier; a two element list containing lower and upper thermal tolerance limits; 
a value indicating trophic level; a value between 0 and 1 indicating trophic specialization; 
the organism's phenotype. Values are sampled from real-world datasets. Trophic level is associated
to specialization, that is we paired the trophic level of a given species from a real food web with
the actual fraction of resources used by that species (over the total number of species in the food
web). In this way, we accounted for non obvious relationships between trophic level and interaction
specialization. Trophic level was computed as the minimum path lenght from target species to basal
resources, thus being an integer number; to that integer part, a random floating part was added
sampled (from a uniform distribution), and used to 'rank' species within the same trophic level
(this will be important to generate realistic food webs, see next section of code). Thermal limits 
were sampled with different frequency from two databases including limits for ectotherm and
endotherm organisms respectively, to account for their different levels of diversity in real systems.
Phenotypes consisted in random strings of letters of size randomly variable in [1,10]. Virtual
'tardigrade' species were also added, with broad thermal tolerance limits and trophic levels
corresponding to those of grazers or micro-predators. Note that not all species generated at this 
stage will be included in the investigated localities/food webs, so the starting global diversity 
in the simulations will be likely lower than the number of species generated here.
 

'''

#generate the set of species
SPN=99900
spp=[]
for i in range(SPN):
	tl,prey_n=sample(tl_file,1)[0]	
	phenotype = ''.join(sample(ft,randrange(1,11)))
	if tl==0:
		ul,ll = sample(plant_file,1)[0]
	else:
		tl+=random()/1000.0
		if random()<0.001:
			ul,ll = sample(endo_file,1)[0]
		else:	
			ul,ll = sample(ecto_file,1)[0]				
	spp.append([i,[ll,ul],tl,prey_n,0,phenotype])	# the 0 value indicate the moment a species arrived in a given locality (the value will change if a species colonizes a locality during the simulation); this is used to compute the time dependent survival probability of a species at a given moment in the simulation)



# add a few tardigrade species
for i in range(100):
	phenotype = ''.join(sample(ft,randrange(1,10)))
	tl,prey_n=sample([j for j in tl_file if j[0] in [1,2]],1)[0]
	tl+=random()/1000.0
	ll,ul = randrange(-80,-20),randrange(50,100)
	spp.append(['tardigrade_'+str(i),[ll,ul],tl,prey_n,0,phenotype])# the 0 value indicate the moment a species arrived in a given locality (the value will change if a species colonizes a locality during the simulation); this is used to compute the time dependent survival probability of a species at a given moment in the simulation)





spp_dict = dict([[i[0],i] for i in spp]) #create a dictionary linking species identifiers to species features





##########################################################################################
### DEFINE FUNCTION TO CREATE FOOD WEBS FROM SPECIES POOL
##########################################################################################


# 1) generate lognormal distribution and rescale it by quoting for maximum value; 
#this will be needed to weigh trophic links in food-webs

'''
2) define food web generating function. Arguments are the species pool (p); the list of possible
functional traits, i.e. alphabet letters (ft); minimum and maximum empirical resource compatibility
(min_tc, max_tc); and the rescaled lognormal distribution (ln_vals). Each entry in the species pool
(p) is a list representing an organism's feature, including: unique identifier; a two element list 
containing lower and upper thermal tolerance limits; a value indicating trophic level; a value between
0 and 1 indicating trophic specialization; the organism's phenotype. Trophic level of a species is
represented by a real number, such as 2.15467. The integer part of that number (2) is the actual trophic
level (those are assigned to species by sampling at random from an empirical distribution of species
trophic levels obtained from real-world food-webs). The floating part of the number (.15467) represent
the relative position of the species within the trophic level, and is used to allow situations such 
as those where two consumers of the same (main) trophic level can prey one on another (with the 
generation of loops in case the two consumers preying one on another share some resource).    
'''

def net_from_pool(p,ft,ft_mat, min_tc, max_tc):
	tl=sorted(list(set([int(i[2]) for i in p])),reverse=True)	# list of trophic levels present in species pool in reverse order (from highest to lowest)
	tl_c=[[i for i in p if int(i[2])==j] for j in tl] #divide species in groups by trophic level
	tl_c_id=[[i[0] for i in j] for j in tl_c] # get the ids of species in each trophic level group
	net=[]
	for i in range(len(tl)-1):	#iterate over trophic levels and create links from one trophic level i to trophic level i+1 (note that the order is reverse, thus links are created from consumers to resources
		for j in tl_c[i]:
			sc=0
			n0=int(round(j[3]*len(p))) #compute the expected number of resources per consumer, on the basis of consumer's specialization, and of the total number of species in the pool
			consumer = j[-1] #j-th consumer's phenotype
			for k in sample(tl_c[i+1],len(tl_c[i+1])): #iterate (in random order) over potential resources of trophic level i+1 (again, this means a trophic level one step lowe than trophic level i)
				if sc<n0:	#check that the maximum number of resources per a given consumer has yet to be reached
					resource = k[-1]	#k-th resource's phenotype
					fc = compute_tc(resource,consumer,ft,ft_mat, min_tc, max_tc) #trait compatibility between consumer and reosurce	
					if random()<fc:
						w = fc	#assign weight to interactions, multiplying trait compatibility per random (rescaled) lognormal value
						net.append([k[0],j[0],w])	#append link to food web
						sc+=1				#count total number of resources used by j-th consumer
			if i<(len(tl)-1):	#check if the expected number of resources (n0) had been assigned to consumer j (on the basis of its specialization)
				for k in sample(tl_c[i],len(tl_c[i])):	#if not replicate above steps, but this time assigning interaction within the same trophic level
					if sc<n0 and j[2]>k[2]:	# as in the description above, comparison is made between the floating parts of the species' trophic level
						resource = k[-1]	#next steps are like those above
						fc = compute_tc(resource,consumer,ft,ft_mat, min_tc, max_tc)	
						if random()<fc:
							w = fc
							net.append([k[0],j[0],w])
							sc+=1
			if i<(len(tl)-2):	#if n0 is yet to be reached, the same procedure as above is repeated, but looking at potential resources distant two trophic levels from the consumer's one.
				for k in sample(tl_c[i+2],len(tl_c[i+2])): #next steps are like those above
					if sc<n0:
						resource = k[-1]
						fc = compute_tc(resource,consumer,ft,ft_mat, min_tc, max_tc)	
						if random()<fc:
							w = fc
							net.append([k[0],j[0],w])
							sc+=1
	g=Graph.TupleList(net,directed=True,edge_attrs="weight") # the food web in form of weighthed edge list is then transformed into a igraph Graph object
	g=g.simplify(combine_edges=max)	#multiple edges are removed, keeping the maximum weight among duplicates
	names=g.vs['name']	#get node names (keeping tracks of names is important because nodes ids in igraph Graphs are different entities than node names; also ids are changed whenever Graph structure is changed
	b_id=set([names.index(i) for i in set(names)&set(tl_c_id[-1])]) #check which nodes in the network are basal resources (trophic level 0)
	to_keep=set([])
	for i in b_id:
		ddd=g.subcomponent(i, mode="OUT") #keep only nodes (species) for which a link to basal resource exists
		to_keep|=set(ddd)
	gg = g.subgraph(to_keep)
	return gg	#return the weighted food-web	



###################################################################################################
###FUNCTIONS TO EVALUATE THE THERMAL COMPATIBILITY BETWEEN A SPECIES AND A LOCALITY
###################################################################################################

#evaluate probability of extinction based on comparison between species environmental tolerance (e.g. thermal) limits, and local conditions (e.g. temperature range)
def niche(sp,loc): #sp and loc are list of two elements; sp includes species thermal tolerance range(minT,maxT);loc includes local thermal range (minT,maxT)  
	if loc[0]<sp[0] or loc[1]>sp[1]: #extinction p is 1 if local T range is not fully contained in species thermal tolerance range
		return 1
	elif min((loc[0]-sp[0]),(sp[1]-loc[1]))>=5:
		return 0
	else:
		return 1/(1+min((loc[0]-sp[0]),(sp[1]-loc[1]))) #otherwise, extinction p remains low for most values of local T, but increases rapidly when local T gets very close to species thermal tolerance limits



def temp_change(t):	#a function describing the relationship between time and temperature 
	return 1*t	#in this case, we assume that one temporal step equals a change of one temperature unit (1 degree)	



def ext_p(sp,loc,t1,t0,scenario,res=100): #compute the probability of a species to be extinct at a given moment; scenario is either 1 for temperature increase, or -1 for temperature decrease
	t = 1+(t1-t0)	#number of temporal units the species has been in a given locality
	loc_0 = loc[0]-(scenario*temp_change(t)), loc[1]-(scenario*temp_change(t)) #identify the environmental conditions (temperature range) of the locality at the beginning of the simulation if target species was already present, or when the species arrived, if that is a colonizer from another locality
	x = linspace(0,t,res)	#identify res = 100 equally spaced moments in between 0 and t
	y = [niche(sp,[loc_0[0]+(scenario*temp_change(x[i])),loc_0[1]+(scenario*temp_change(x[i]))]) for i in range(res)] #identify target species' survival probability at each moment (based on species' niche and environmental conditions at a given moment).
	return trapz(y,x)	#compute the cumulative probability of extinction as the area under the curve y=f(x)



###################################################################################################
### DEFINE FUNCTION TO SIMULATE CO-EXTINCTION CASCADES
###################################################################################################

'''
The function takes as input (g): a food web in form of an igraph Graph object. (eee): an array including
a list of species ids going extinct (primary extinctions). (spp_dict): a dictionary mapping species 
ids to species features; (pred_0_dict): a dictionary mapping species ids to the starting amount of
 resources used by a consumer,i.e. the cumulated weight of directed resource-consumer interactions 
for the target consumer; those include all trophic links where the consumer acts as a consumer, 
not those where it acts as a resource; starting amount means from the beginning of time, not from the
beginning of the current extinction cascade (prey_0_dict): a dictionary mapping species ids to the 
starting amount of consumers using a certain resource, i.e. the sum of all directed resource-consumer
interactions involving target resource. (tre): a threshold value indicating the minimum fraction of
lost resources triggering the secondary extinction of a consumer; such fraction is measured, at each
step of the co-extinction cascade, as the ratio between the starting total amount of consumer-resource
interactions involving the target consumer (of which track is kept by pred_0_dict) and the reamaining
amount of interactions at the given step.(rew): the fraction of resources made available by the 
extinction of a consumer that can be reallocated to other consumers.     
'''


def ext_casc(g,eee,pred_0_dict,prey_0_dict,spp_dict,tre=0.5,rew=0.5):
	names = g.vs['name']
	gw_mat = g.get_adjacency(attribute='weight')	#igraph Graph object is transformed in a weighted adjacency matrix 
	gw_mat = array([[el for el in row] for row in gw_mat],dtype = float) #matrix is converted to array
	eee = [names.index(i) for i in eee]	#get the position in the array of the nodes going extinct 
	new_mat=gw_mat.copy()	#create a copy of the original matrix
	pp = set([i for i in range(len(names)) if spp_dict[names[i]][2] == 0])	#identify basal resources, as species having trophic level 0
	pred_0 = array([pred_0_dict[i] for i in names]) #starting amount of resources using by consumers in the food web
	prey_0 = array([prey_0_dict[i] for i in names])	#starting amount of consumers using the resources in the food web
	ext=set(eee)	#transform array to set
	col_sum_0 = gw_mat.sum(0) #initial column totals; those correspond to the amount of resources used by consumers before the  co-extinction cascade starts	
	while len(eee)>0:	#reiterate the steps in the co-extinction cascade as long as new secondary extinctions are triggered
		new_mat[eee,:]=0.0	#primary extinctions step 1: delete interactions where the extinct species act as consumers 
		new_mat[:,eee]=0.0	#primary extinctions step 2: delete interactions where the extinct species act as resources
		col_sum_diff  = new_mat.sum(0)/col_sum_0	#compute ratio between resources used by consumer after and before the co-extinction step		
		col_sum_diff[where(col_sum_0==0)] = 1.0	#set to 1 the values in the ratio for which the initial value of used resources was 0
		new_mat = (col_sum_diff*new_mat.T).T	#the reduction in available resources for a given consumer is projected to the consumer's consumers, under the assumption that a reduction in available resources used by a consumer would result in a reduction of consumer's population
		eee=set(where(new_mat.sum(0)/pred_0<=tre)[0])-(ext)	#identify secondary extinctions, as consumer losing a fraction of resources equal or larger than tre; this is evaluated on the basis of the resource used at the beginning of time, not at the beginning of the co-extinction cascade
		ext|=eee	#kept track of cumulated extinctions
		eee=array(list(eee))
		prey_avail=prey_0-new_mat.sum(1) #evaluate the amount of resources that have been freed by primary extinctions
		prey_avail[array(list(ext))]=0.0 #exclude from those the resources that have gone extinct
		pred_n = (new_mat>0).sum(1)	#compute the numbers of consumers using each resource in the network 
		to_add=((prey_avail/pred_n)*rew)	#compute the amount of resources that can be reallocated to consumers; note that resources are reallocated only to species that are already using them, i.e. the structure of the network is not changed, but the weights of the consumers receiving the reallocable resource/s increase; allocable resource is equally divided between recipient consumers, and reduced by a 'rewiring' factor (rew).
		to_add[where(pred_n==0)]=0	#replace na values with 0
		to_add=(to_add*(new_mat.T>0))	#filter the to_add array to ensure the resource will be added to the proper set of consumers (i.e. those already using the resource, see previous comments)
		new_mat=(new_mat.T+to_add).T	#reallocate the resource
	g = Graph.Weighted_Adjacency([list(row) for row in gw_mat])	#convert the adjacency matrix back to an igraph Graph object
	g.vs['name'] = names	
	g = g.subgraph(set(range(len(g.vs)))-ext)	#in the conversion process, nodes with 0 degree (i.e. with no in- or out-going links) are included in the network; this step takes them out
	ext_names = set([names[i] for i in ext]) #get the names (i.e. unique ids) of all extinct nodes/species
	b_id=[node.index for node in g.vs if spp_dict[node['name']][2]==0] #identify basal resources (i.e. species in the food web having trophic level = 0)
	pre_check = set(g.vs['name'])	#get names of all nodes in the new network before the final step, where all nodes having no path to basal resources are deleted from the network
	to_keep=set([])
	for i in b_id:
		ddd=g.subcomponent(i, mode="OUT") #keep only species for which a path to basal resource exists
		to_keep|=set(ddd)
	g = g.subgraph(to_keep)	#reduce the graph to only species having a path connecting them to basal resources
	post_check = set(g.vs['name'])	#get the names of the reduced graph
	ext_names|=(pre_check-post_check) #add the species that have been removed by this last step to the set of extinct species
	return ext_names, g	#return the full set of species going extinct following the cascade, and the resulting food-web 




#################################################################################################### 
###GENERATE LOCALITIES, POPULATE THEM WITH SPECIES, BUILD FOOD WEBS
####################################################################################################

pools = []	#empty list of species pools; each pool is the list of species found in a given locality
nets = []	#empty list of food webs
space = []	#empty list of localities' position and climate (lat,lon,minT,maxT)
glob_div = randrange(100,1000)	#number of species sampled at random in the global diversity set to be 'dropped' in each locality
loc_N = randrange(100,500)
while len(space)<loc_N:	#number of localities to be populated with species
	x,y,mean_T,min_T,max_T=sample(env_file,1)[0]	#sample lat,lon,minT,maxT to be assigned to new locality; those are sampled from real climatic data
	if [x,y] not in [j[0] for j in space]:	#check that the new locality does not overlap with existing ones
		loc=[[x,y],[min_T,max_T],'na']	#assign the sampled features to the new locality
		pool=[]		#generate empty pool	
		for attempt in range(glob_div):	#drop the selected number of random species in the new locality, one after another
			sp=sample(spp,1)[0]	#sample a candidate species
			if sp not in pool and random()>niche(sp[1],loc[1]): #check compatibility between the species thermal limits and local minT and maxT
				pool.append(sp)	#if the species is compatible, add it to the pool
		loc[2]=len(pool)	#NOT NEEDED?	
		if len(pool)>=5 and len([i for i in pool if i[2]==0])/float(len(pool))>=0.20:	#if enough species are in pool (>4) and there are enough basal resources (>=20% of local diversity), try to build a food web
			net = net_from_pool(pool,ft,ft_mat, min_tc, max_tc)
			if len(net.es)>0: 	#if the attempt to build a food web succeeds:
				nets.append(net)	#add the web to the list of networks;
				in_net=set(net.vs['name'])	#check which species are in the network
				pool = [j for j in pool if j[0] in in_net]	#include in pool only species that have links in the network
				pools.append(pool)	#add local pool to pools list
				space.append(loc)	#add the new locality to the locality list
				print len(space),len(pool),loc[1],len([i for i in pool if i[2]==0])/float(len(pool)) #print some data




#############################################################################################################################
####GENERATE DICTIONARIES MAPPING THE INITIAL AMOUNT OF RESOURCES USED BY A CONSUMER, AND OF CONSUMERS USING A GIVEN RESOURCE
#############################################################################################################################

'''
Amount of resources used by consumers, and of consumers using resources is quantified in terms of 
cumulative weight of interactions. The dictionaries are needed by the co-extinction cascade function
(ext_casc). Note that one dictionary is generated for each food-web.
'''

preds_0, preys_0 = [],[]
for net in nets:	#iterate over all food-webs
	wg = net.get_adjacency(attribute='weight')	#convert food web to weighted adjacency matrix
	wg = array([[cell for cell in row] for row in wg],dtype = float)	#convert matrix to array
	names = net.vs['name']	#get names from network
	spn = len(names)
	pred_sum = wg.sum(0) #get sums of columns, i.e. total amounts of resources used by consumers	
	prey_sum = wg.sum(1)	#get sums of rows, i.e. total amounts of consumers using a resource
	preds_0.append(dict([[names[i],pred_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	preys_0.append(dict([[names[i],prey_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	
	





####################################################################################################
###PRELIMINARY DISPERSAL/COLONIZATION PHASE
####################################################################################################

'''
This is meant to add a bit of biogeographical meaning to the communities; species are moved at 
random between localities, with colonization happening with probability determined by distance
between source and target locality, network structure (edge density) at the target locality,
as well as the adaptation of the colonizer to the environmental conditions (i.e. temperature range)
of the target locality, and its ability to replace (or co-exist with) species of the food web of 
the target locality.
'''
eta=-1	#exponent of the dispersal kernel 
kkk=1	#constant value of the dispersal kernel
pre_col_att_n = randrange(1000,100000) #set the number of colonization attempts
for col_att in range(pre_col_att_n):	
	i,j=sample(range(len(pools)),2)	#sample target and source localities 
	if len(nets[i].es)*len(nets[j].es)>0:	#check that both localities host some species embedded in a food web
		invasability = 1-len(nets[j].es)/float(len(nets[j].vs)**2) #determine network invasability as 1 minus edge density (given by the total number of links divided by the squared number of nodes)
		d=euclidean(space[i][0],space[j][0])	#measure distance between source and target locality
		pd=kkk*d**eta	#evaluate probability of dispersal from source to target locality, decreasing exponentially with distance
		if random()<pd:	#colonization attempts take place only with probability pd
			col_sp=sample(pools[i],1)[0]	#sample a potential colonizer from locality i
			names = nets[j].vs['name']	#get names of food web of locality j
			new_list = []	#create empty list for new links
			col_ok = 'no'	#label for successful colonization
			for tlink in nets[j].es:	# iterate over trophi links in food web of locality j
				edge_replaced = 'no'	#label for successful edge (i.e. trophic link) replacement						
				prey,pred = tlink.tuple	#first node of an edge is resource, second is consumer
				if random()<=invasability:	#colonization is constrained by network invasability (see above)
					if int(spp_dict[pred][2]) == int(col_sp[2]):	#evaluate if the colonizer can replace a consumer in using a certain resource; check if the trophic level of the colonizer is the same as that of the consumer in the edge/link under evaluation 
						if niche(col_sp[1],space[j][1])<=niche(spp_dict[pred][1],space[j][1]): #colonization can happen only if the colonizer has equal or higher compatibility to local environmental conditions than the consumer it attempts to replace
							if random()<=compute_tc(spp_dict[prey][-1],col_sp[-1],ft,ft_mat, min_tc, max_tc): #a last constraint is given by the trait compatibility between the colonizer and the target resource; for this, the two phenotypes are compared
								new_list.append([names[prey],col_sp[0],tlink['weight']]) #if all the previous requirements are met, colonization succeeds, and a link is added to the new edge list ('new_list'), connecting the original resource to the colonizer, with the original weight of the interaction 
								edge_replaced = 'yes'	#change status of successful edge replacement label
								col_ok = 'yes'	#change status of successful colonization label
					elif int(spp_dict[prey][2]) == int(col_sp[2]):	#evaluate if the colonizer can become a resource to other species; check if the trophic level of a resource is the same as that of the colonizer
						if random()<=compute_tc(col_sp[-1],spp_dict[pred][-1],ft,ft_mat, min_tc, max_tc): #check trait compatibility
							if int(col_sp[2])==0 and niche(col_sp[1],space[j][1])<=niche(spp_dict[prey][1],space[j][1]):	#if the colonizer is a basal resource better adapted to local climate than the local resource, it will outcompete it
								new_list.append([col_sp[0],names[pred],tlink['weight']]) #if all the previous requirements are met, colonization succeeds, and a link is added to the new edge list ('new_list'), connecting the original consumer to the colonizer, with the original weight of the interaction 
								edge_replaced = 'yes' #change status of successful edge replacement label
								col_ok = 'yes' #change status of successful colonization label
							elif int(col_sp[2])>0 and random()>niche(col_sp[1],space[j][1]): #if the colonizer is not a basal resource and adapted to the climate of the target locality, it will be used as an alternative resource by the consumer (not outcompeting local resource)
								new_list.append([names[prey],names[pred],tlink['weight']*random()]) #add the new link between the colonizer as a resource and the local consumer; weight of the new interaction is a random fraction of the weight of the target interaction; this actually strenghten the food web
								new_list.append([col_sp[0],names[pred],tlink['weight']])	#add the pre-existing link, with the original weight
								edge_replaced = 'yes'	#change status of successful edge replacement label
								col_ok = 'yes'	#change status of successful colonization label
				if edge_replaced!= 'yes':	#if no colonization attempt was successful, add back the original link to the network
					new_list.append([names[prey],names[pred],tlink['weight']])
			new_net = Graph.TupleList(new_list,directed=True,edge_attrs='weight')# after all network links have been evaluated for potential colonization, convert the edge list into an igraph Graph object
			new_net.simplify(combine_edges=sum)	#multiple edges are collapsed into single ones, and the weights of multiple edges are summed up 
			if col_ok != 'no':	#in case of successful colonization, the dictionaries listing the initial amount of resources used by a consumer, and of consumers using a resource need to be updated
				new_id = new_net.vs['name'].index(col_sp[0]) #get the node index of the colonizer in the new network
				prey_0_new = new_net.strength(new_id,weights='weight',mode='OUT') #get the total amount of resources used by the colonizer
				pred_0_new = new_net.strength(new_id,weights='weight',mode='IN') #get the total amount of consumers using the colonizer	
				preys_0[j][col_sp[0]] = prey_0_new	#update the dictionary of the initial amount of resources used by consumers; note that the initial amount for the colonizer is considered here from the moment of the successful colonization 
				preds_0[j][col_sp[0]] = pred_0_new	#update the dictionary of the initial amount of consumers using each resource; note that the initial amount for the colonizer is considered here from the moment of the successful colonization
				ext_,new_net = ext_casc(new_net,[],preds_0[j],preys_0[j],spp_dict,tre=0.5,rew=0.5) #perform an extinction  cascade with an empty set of primary extinction; this is made as a shortcut to prune the network from nodes not having links to basal resources
				nets[j] = new_net #update the network in locality j
				pools[j] = [spp_dict[new_net.vs['name'][k]] for k in range(len(new_net.vs))] #update species pool in locality j





#############################################################################################################################
####RESET DICTIONARIES MAPPING THE INITIAL AMOUNT OF RESOURCES USED BY A CONSUMER, AND OF CONSUMERS USING A GIVEN RESOURCE
#############################################################################################################################

'''
This step reset the two dictionaries mapping the amount of resources used by consumers, and of consumers using resources is quantified in terms of 
cumulative weight of interactions. Update is needed because we assume that the networks are at equilibrium before the environmental change
phase.
'''

preds_0, preys_0 = [],[]
for net in nets:	#iterate over all food-webs
	wg = net.get_adjacency(attribute='weight')	#convert food web to weighted adjacency matrix
	wg = array([[cell for cell in row] for row in wg],dtype = float)	#convert matrix to array
	names = net.vs['name']	#get names from network
	spn = len(names)
	pred_sum = wg.sum(0) #get sums of columns, i.e. total amounts of resources used by consumers	
	prey_sum = wg.sum(1)	#get sums of rows, i.e. total amounts of consumers using a resource
	preds_0.append(dict([[names[i],pred_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	preys_0.append(dict([[names[i],prey_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	
	


###########################################################################################################################
###PERFORM GLOBAL CHANGE EXPERIMENT TO ASSESS THE RELATIVE IMPORTANCE OF CO-EXTINCTIONS IN GLOBAL BIODIVERSITY DECLINE
###########################################################################################################################
'''
Global change is simulated as a widespread, monotonic increase (or decrease) in local temperature maxima (or minima).
The change in temperature triggers the local extinction of species no longer able to survive in the novel conditions.
In a first scenario, such primary extinctions are considered not able to trigger secondary extinctions; in a second 
scenario, primary extinctions are permitted to trigger co-extinction cascades. In the comments, the first scenario is 
referred to as the 'tolerance' scenario, while the second one as the 'co-extinction' scenario. The Diversity decline 
is tracked in both scenarios, over 100 steps of temperature change. In between climate change steps, movements of 
species between localities is permitted in the same way as in the preliminary colonization phase. In the scenario not 
accounting for co-extinctions, the process is simplified by not accounting for network structure in determining the 
success of a colonization event. Species are also allowed to (partially) adapt to changing conditions, by progressively 
increasing their tolerance to thermal extremes. Furthermore, a rescue effect is implemented, as the stochastic arrival 
of novel recruits from neighboring localities, replenishing depleted resources. 
'''

all_spp = set([item for sublist in [[i[0] for i in j] for j in pools]  for item in sublist]) #get global species pool 
biodiv_0 =  float(len(all_spp))#measure overall starting biodiversity
biodiv_0_tard = float(len([i for i in all_spp if 'tardigrade' in str(i)])) #measure overall starting tardigrade biodiversity

 
scenario = sample([1,-1],1)[0] # this value indicate the direction of the temperature change; value 1 indicates heating, set the value to -1 for planetary cooling 
tre_val = random() #set the co-extinction treshold value (the 'tre' parameter in the 'ext_casc' function); this could be set, for example to [0.25,0.50,0.75]
rew_val = random() #set the resource reallocation ratio (the 'rew' parameter in the 'ext_casc' function); this could be set, for example to [0.25,0.50,0.75]
adpt_p = random()/100.0	#determine the individual probability of a species to adapt its thermal niche (i.e. either increase its upper temperature limit, or decrease the lower one)
col_att_n = randrange(100,10000)
rescue_prob = 1

fname = str(random())[2:]+'_'+'_'.join(map(str,[biodiv_0,len(space),pre_col_att_n,col_att_n,tre_val,rew_val,adpt_p,scenario]))+'.csv'

pools_1=deepcopy(pools)	#create a copy of the pools list for the 
pools_2=deepcopy(pools)
nets_1=deepcopy(nets)	#igraph weighted networks
new_space=deepcopy(space)	#for the rescue effect
RES=[[0.0,1.0,1.0,1.0,1.0]]
out=open('./results/'+fname,'w')	#open output file in append mode
out.write(','.join(map(str,RES[-1]))+'\n')	#write results to file
out.close()	#close output file


###ACCOUNTING FOR LATITUDINAL EFFECT
#we assume that climate change could happen faster at higher latitude (>60 degrees N/S); based on ongoing change and current knowledge (e.g. https://www.nature.com/articles/ngeo2234), we assumed a linear increase between 60 and 90 degrees, with probabilit of temperature change double at poles
def lat_effect(lat):
	if lat<60:
		lat_factor = 1
	else:
		lat_factor = 1+ (abs(lat)-60)/30.0
	return  lat_factor #this is a factor to be applied to loc change 



T = 0 # temperature change
time_res = 100.0
for time_step in range(int(50*time_res)): #set number of step
	av_change = []	#list to be filled with local T changes, needed to compute global average change
	exts_1,exts_2=[],[]	#empty lists to be filled with the sets of extinct species ids in the tolerance and co-extinction scenarios in each locality				
	for i in range(len(pools)): #iterate over localities/communities
		change_val = abs(normal(1.0,0.25))/time_res	#assign a change in local temperature sampled from a normal distribution
		lat_factor = lat_effect(new_space[i][0][0])
		av_change.append(change_val)	#keep track of local change
		new_space[i][1][0]+=(scenario*change_val*lat_factor) #change min T 
		new_space[i][1][1]+=(scenario*change_val*lat_factor) #change max T
		ext_1,ext_2=[],[]	#identify primary extinctions due to T change
		for j in pools_1[i]:	# evaluate the compatibility between all species in the pool with the new climate in the tolerance scenario
			if random()<1/time_res:	
				if random()<ext_p(j[1],new_space[i][1],time_step/time_res,j[-2],scenario,res=100):	
					ext_1.append(j[0])
		for j in pools_2[i]:	# do the same for the co-extinction scenario; evaluation is made separately, because species pools will clearly become different one from another throughout the simulation
			if random()<1/time_res:
				if random()<ext_p(j[1],new_space[i][1],time_step/time_res,j[-2],scenario,res=100):
					ext_2.append(j[0])
		exts_1.append(list(set(ext_1)))	#append the set of extinct species ids in the tolerance scenario to the global list (including sets of extinct species for all localities)
		exts_2.append(list(set(ext_2)))	# do the same for the co-extinction scenario
	alive_1=set([])	#create an empty set to contain the ids of all extant species in the tolerance scenario
	for i in range(len(pools)):	#count the number of extant species in the tolerance scenario; for this, iterate through the pools 
		pools_1[i]=[j for j in pools_1[i] if j[0] not in exts_1[i]] #eliminate species that went extinct (primary extinctions)
		alive_1|=set([j[0] for j in pools_1[i]])	#add the other species ids to the global set of extant species 
	tard_1=len(set([j for j in alive_1 if 'tardigrade' in str(j)]))	#count the number of alive tardigrade species
	alive_2=set([])	#create an empty set to contain the ids of all extant species in the co-extinction scenario	
	for i in range(len(pools)): #iterate over all localities
		pools_1[i]=[j for j in pools_1[i] if j[0] not in exts_1[i]]	#eliminate extinct species from local pool in the tolerance scenario
		pools_2[i]=[j for j in pools_2[i] if j[0] not in exts_2[i]]	# do the same for the co-extinction scenario
		if len(nets_1[i].es)>0:	#check if local food web has extant links
			names=nets_1[i].vs['name']
			ext_3,new_g=ext_casc(nets_1[i],exts_2[i],preds_0[i],preys_0[i],spp_dict,tre=tre_val,rew=rew_val) #evaluate co-extinction cascades triggered by primary extinctions in locality i
			pools_2[i]=[j for j in pools_2[i] if j[0] not in ext_3] #reduce species pool to species still in the network after the co-extinction cascade
			nets_1[i]=new_g #replace the starting network with the one obtained after the co-extinction cascade
		alive_2|=set([k[0] for k in pools_2[i]]) #add the set of alive species ids to the global one
	tard_2=len(set([j for j in alive_2 if 'tardigrade' in str(j)]))	#count the number of alive tardigrades (in the co-extinction scenario)
	T+=mean(av_change) #keep track of the average increase in T
	###rescue mechanism
	for i in range(len(pools)): #iterate over all pairwise combinations of different localities
		for j in range(len(pools)):
			if i!=j and random()<rescue_prob/time_res:
				d=euclidean(space[i][0],space[j][0]) #compute distance between localities
				pd=kkk*d**eta	#compute distance based probability of movement from source to target locality
				pool_j = [sp[0] for sp in pools_2[j]]	#list of species in locality j
				names_i = nets_1[i].vs['name']	#names of species in network i (and hence in locality i)
				pred_gap = [preds_0[i][names_i[k]]-nets_1[i].strength(k,weights='weight',mode='IN') for k in range(len(nets_1[i].vs))] #compute the amount of each resource that have been depleted for a consumer from the beginning of the simulation; this is, for each consumer, the difference beween the original cumulative weight of all interactions involving the target consumer, and the same quantity at the given simulation step 
				new_list = [] #create an empty list of links
				time_reset = [] #list of species for which time in a locality must be reset
				for tlink in nets_1[i].es:	#iterate over all links in the network
					prey,pred = tlink.tuple	#each link is composed by a resource and its consumer
					if random()<=pd and names_i[prey] in pool_j:	#the resource in the target link is replenished with a probability pd determined by distance between source and target locality, and only if the resource in the target link exists also in the source locality (j)
						add_weight = tlink['weight']*pd #if the conditions in the previous line are met, the weight of the target link is increased by its own value multiplied by the pd, under the assumption that the exchange of novel recruits will be stronger for closer localities
						if add_weight>pred_gap[pred]: #the rescuing is thought just as a mechanism to replenish depleted resources for a consumer; thus the cumulative interaction weight for a consumer cannot exceed the starting one (i.e. that before climate change); this step takes care of this 
							add_weight = pred_gap[pred] #the maximum interaction weight that can be added is given by the difference between the cumulative original interaction weight of target consumer (i.e. before climate change), and the current one (taking into account other potential rescue events happened in this same simulation step). 
						pred_gap[pred]-=add_weight	#update the amount of depleted resource for the target consumer
						new_list.append([names_i[prey],names_i[pred],add_weight+tlink['weight']]) #add the replenished interaction to the edge list
						time_reset.append(names_i[prey])	#list of local species that have received new recruits; for those, the time of exposition to local environmental conditions (used to compute survival probability at a given moment) is reset   
					else:
						new_list.append([names_i[prey],names_i[pred],tlink['weight']]) # with probability 1-pd and/or if the target resource is not in species pool of locality j, the original edge is added to the edge list with no changes
				new_net = Graph.TupleList(new_list,directed=True,edge_attrs='weight') #rebuild the igraph Graph network from edgelist
				new_net.add_vertices(list(set(names_i)-set(new_net.vs['name'])))	#add back to the network nodes with 0 degree (those will be removed by the next co-extinction cascade, but need to be there to avoid errors in the code execution)
				nets_1[i] = new_net #update network i
				for k in range(len(pools_2[i])):	#iterate over species in the target pool, and reset time of permanence/exposition to local environmental conditions for those species that have received new recruits
					if pools_2[i][k][0] in time_reset:
						pools_2[i][k][-2] = (time_step)/time_res
	###dispersal/colonization phase in the co-extinction scenario (same as in the preliminary dispersal/colonization phase)
	for col_att in range(col_att_n):
		i,j=sample(range(len(pools)),2)
		#determine invasability of net j
		if len(nets_1[i].es)*len(nets_1[j].es)>0:
			invasability = 1-len(nets_1[j].es)/float(len(nets_1[j].vs)**2) #invasability depends on network connectance
			d=euclidean(space[i][0],space[j][0])	#measure distance between source and target locality
			pd=kkk*d**eta	#compute dispersal probability based on distance
			if random()<pd/time_res:	#colonization is attempted on the basis of pd; the number of colonization attempts per step depends also on the temporal resolution of the simulation (i.e. the time interval between subsequent steps)
				col_sp=sample(pools_2[i],1)[0][:]	#sample a potential colonizer from locality i; note that a copy of the colonizer is created (by adding '[:]' at the end of the line); this is important to avoid that time of permanence is reset for the colonizer in both the source and the target locality (it needs to be reset only in the latter)
				if col_sp[0] not in [k[0] for k in pools_2[j]]: #the arrival of new recruits of a species already present are already accounted for by the rescue mechanism
					names = nets_1[j].vs['name']	#get names of food web of locality j
					new_list = []	#create empty list for new links
					col_ok = 'no'	#label for successful colonization
					for tlink in nets_1[j].es:	# iterate over trophi links in food web of locality j
						edge_replaced = 'no'	#label for successful edge (i.e. trophic link) replacement						
						prey,pred = tlink.tuple	#first node of an edge is resource, second is consumer
						if random()<=invasability:	#colonization is constrained by network invasability (see above)
							if int(spp_dict[pred][2]) == int(col_sp[2]):	#evaluate if the colonizer can replace a consumer in using a certain resource; check if the trophic level of the colonizer is the same as that of the consumer in the edge/link under evaluation 
								if niche(col_sp[1],space[j][1])<=niche(spp_dict[pred][1],space[j][1]): #colonization can happen only if the colonizer has equal or higher compatibility to local environmental conditions than the consumer it attempts to replace
									if random()<=compute_tc(spp_dict[prey][-1],col_sp[-1],ft,ft_mat, min_tc, max_tc): #a last constraint is given by the trait compatibility between the colonizer and the target resource; for this, the two phenotypes are compared
										new_list.append([names[prey],col_sp[0],tlink['weight']]) #if all the previous requirements are met, colonization succeeds, and a link is added to the new edge list ('new_list'), connecting the original resource to the colonizer, with the original weight of the interaction 
										edge_replaced = 'yes'	#change status of successful edge replacement label
										col_ok = 'yes'	#change status of successful colonization label
							elif int(spp_dict[prey][2]) == int(col_sp[2]):	#evaluate if the colonizer can become a resource to other species; check if the trophic level of a resource is the same as that of the colonizer
								if random()<=compute_tc(col_sp[-1],spp_dict[pred][-1],ft,ft_mat, min_tc, max_tc): #check trait compatibility
									if int(col_sp[2])==0 and niche(col_sp[1],space[j][1])<=niche(spp_dict[prey][1],space[j][1]):	#if the colonizer is a basal resource better adapted to local climate than the local resource, it will outcompete it
										new_list.append([col_sp[0],names[pred],tlink['weight']]) #if all the previous requirements are met, colonization succeeds, and a link is added to the new edge list ('new_list'), connecting the original consumer to the colonizer, with the original weight of the interaction 
										edge_replaced = 'yes' #change status of successful edge replacement label
										col_ok = 'yes' #change status of successful colonization label
									elif int(col_sp[2])>0 and random()>niche(col_sp[1],space[j][1]): #if the colonizer is not a basal resource and adapted to the climate of the target locality, it will be used as an alternative resource by the consumer (not outcompeting local resource)
										new_list.append([names[prey],names[pred],tlink['weight']*random()]) #add the new link between the colonizer as a resource and the local consumer; weight of the new interaction is a random fraction of the weight of the target interaction; this actually strenghten the food web
										new_list.append([col_sp[0],names[pred],tlink['weight']])	#add the pre-existing link, with the original weight
										edge_replaced = 'yes'	#change status of successful edge replacement label
										col_ok = 'yes'	#change status of successful colonization label
						if edge_replaced!= 'yes':	#if no colonization attempt was successful, add back the original link to the network
							new_list.append([names[prey],names[pred],tlink['weight']])
					new_net = Graph.TupleList(new_list,directed=True,edge_attrs='weight')# after all network links have been evaluated for potential colonization, convert the edge list into an igraph Graph object
					new_net.simplify(combine_edges=sum)	#multiple edges are collapsed into single ones, and the weights of multiple edges are summed up 
					if col_ok != 'no':	#in case of successful colonization, the dictionaries listing the initial amount of resources used by a consumer, and of consumers using a resource need to be updated
						new_id = new_net.vs['name'].index(col_sp[0]) #get the node index of the colonizer in the new network
						prey_0_new = new_net.strength(new_id,weights='weight',mode='OUT') #get the total amount of resources used by the colonizer
						pred_0_new = new_net.strength(new_id,weights='weight',mode='IN') #get the total amount of consumers using the colonizer	
						preys_0[j][col_sp[0]] = prey_0_new	#update the dictionary of the initial amount of resources used by consumers; note that the initial amount for the colonizer is considered here from the moment of the successful colonization 
						preds_0[j][col_sp[0]] = pred_0_new	#update the dictionary of the initial amount of consumers using each resource; note that the initial amount for the colonizer is considered here from the moment of the successful colonization
						ext_,new_net = ext_casc(new_net,[],preds_0[j],preys_0[j],spp_dict,tre=0.5,rew=0.5) #perform an extinction  cascade with an empty set of primary extinction; this is made as a shortcut to prune the network from nodes not having links to basal resources
						nets_1[j] = new_net #update the network in locality j					
						pools_2[j] = [k for k in pools_2[j] if k[0] in new_net.vs['name']] #update species pool in locality j, by excluding species not in network
						if col_sp[0] in new_net.vs['name']:	
							col_sp[-2] = time_step/time_res		#reset time of permanence/exposition to local conditions for the colonizer							
							pools_2[j].append(col_sp)	#add the successful colonizer to species pool
		###dispersal/colonization in the tolerance scenario
		i,j=sample(range(100),2)	#select two localities at random
		if len(pools_1[i])>0:	#check that the source locality has some species
			invasability = 1-(len(pools_1[i])/float(len(pools[i]))) #invasability increases with diversity loss
			d=euclidean(space[i][0],space[j][0])	#measure distance between source and target locality
			pd=kkk*d**eta	#compute probability of dispersal based on distance betweeen localities
			if random()<pd/time_res:	#colonization is attempted with probability pd
				col_sp=sample(pools_1[i],1)[0][:] #select at random a potential colonizer; note that a copy of the colonizer is created (by adding '[:]' at the end of the line); this is important to avoid that time of permanence is reset for the colonizer in both the source and the target locality (it needs to be reset only in the latter)
				if col_sp[0] not in [k[0] for k in pools_1[j]]: #the arrival of new recruits of a species already present are already accounted for by the rescue mechanism
					to_del = []	#create empty list to be populated with species in target locality out-competed by the colonizer
					col_ok = 'no' #label for successful colonization
					competitors = 'no' #label indicating the presence of competitors of the colonizer in target locality							
					for loc_sp in pools_1[j]:	#iterate over species in target locality 
						if int(loc_sp[2]) == int(col_sp[2]) and random()<phen_similarity(col_sp[-1],loc_sp[-1]): #evaluate the similarity between the functional traits of the colonizer and that of the local species, if both belong to the same trophic level
							competitors = 'yes' #change status of competitors label	
							if random()<invasability and niche(col_sp[1],new_space[j][1])<=niche(loc_sp[1],new_space[j][1]):#the colonizer outcompete a local species of the same trophic level and with similar functional traits if better adapted than the latter to local climate and with probability determined by invasability of target locality
								to_del.append(loc_sp[0]) #if above conditions are met, colonization succeeds
								col_ok = 'yes' #change status of successful colonization label
					if col_ok == 'no' and competitors =='no': #colonization can happen without outcompetition if colonizer can enter an unoccupied niche
						if random()<invasability and random()>niche(col_sp[1],new_space[j][1]): #in this case, colonization success is determined by local invasability, and by compatibility between colonizer's niche and local climate 
							col_ok = 'yes'	#change status of successful colonization label
					if to_del != []:	#if any local species has been outcompeted by the colonizer, species pool in j need to be updated
						pools_1[j] = [sp for sp in pools_1[j] if sp[0] not in to_del] #outcompeted species are eliminated from the j-th pool
					if col_ok != 'no':								
						col_sp[-2] = time_step/time_res	#reset time of permanence/exposition to local conditions for the colonizer
						pools_1[j].append(col_sp)	#the colonizer is added to the j-th pool
	### dynamical species adaptation to changing conditions
	for pool in range(len(pools)):	#iterate over all localities/communities
		for sp in range(len(pools_1[pool])):	#iterate over all species in a locality in the tolerance scenario
			if random()<adpt_p/time_res:	#adaptation happens with probability adpt_p
				adpt = scenario*abs(normal(0.75,0.25))	#the extent of thermal adaptation is a random value sampled from a normal distribution; sign depends on the scenario (heating or cooling)
				pools_1[pool][sp][1][0]+=adpt	#the change in tolerance is applied to the lower species thermal limit
				pools_1[pool][sp][1][1]+=adpt	#the change in tolerance is applied to the upper species thermal limit
		for sp in range(len(pools_2[pool])):	#the same steps are applied to the species pools in the co-extinction scenario
			if random()<adpt_p/time_res:
				adpt = scenario*abs(normal(0.75,0.25))
				pools_2[pool][sp][1][0]+=adpt
				pools_2[pool][sp][1][1]+=adpt
	#keep track of results, and write them to file
	RES.append([T,len(alive_1)/biodiv_0,len(alive_2)/biodiv_0,tard_1/biodiv_0_tard,tard_2/biodiv_0_tard]) #tracked results include average T change, totalt extant diversity in the tolerance scenario, total extant diversity in the co-extinction scenario, and the numbers of extant 'tardigrade' species in both scenarios
	print ','.join(map(str,RES[-1]))#print results on screen
	out=open('./results/'+fname,'a')	#open output file in append mode
	out.write(','.join(map(str,RES[-1]))+'\n')	#write results to file
	out.close()	#close output file








#############################################################################################################################
####EVALUATE NETWORK ROBUSTNESS TO WARMING AND COOLING
#############################################################################################################################

'''
First, reset the two dictionaries mapping the amount of resources used by consumers, and of consumers using resources is quantified in terms of 
cumulative weight of interactions. Update is needed because we assume that the networks are at equilibrium before the environmental change
phase.
'''

preds_0, preys_0 = [],[]
for net in nets:	#iterate over all food-webs
	wg = net.get_adjacency(attribute='weight')	#convert food web to weighted adjacency matrix
	wg = array([[cell for cell in row] for row in wg],dtype = float)	#convert matrix to array
	names = net.vs['name']	#get names from network
	spn = len(names)
	pred_sum = wg.sum(0) #get sums of columns, i.e. total amounts of resources used by consumers	
	prey_sum = wg.sum(1)	#get sums of rows, i.e. total amounts of consumers using a resource
	preds_0.append(dict([[names[i],pred_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	preys_0.append(dict([[names[i],prey_sum[i]] for i in range(spn)]))	#make dictionary mapping values to species ids
	
	


	
pools_1=deepcopy(pools)	#create a copy of the pools list for the 
nets_1=deepcopy(nets)	#igraph weighted networks


out=open('./results/'+'rob_'+fname,'w')	#open output file in append mode


for p in range(len(pools)):
	ext_seq_warm = [j[2] for j in sorted([[i[1][1],random(),i[0]] for i in pools_1[p]])]
	ext_seq_cool = [j[2] for j in sorted([[i[1][0],random(),i[0]] for i in pools_1[p]])[::-1]]
	n = float(len(ext_seq_warm))
	node_deg = [[i.degree(mode='OUT'),random(),i['name']] for i in nets_1[p].vs]
	ext_seq_best = [i[2] for i in sorted(node_deg)]
	ext_seq_worst = [i[2] for i in sorted(node_deg)[::-1]]
	ext_seq_random = sample(nets_1[p].vs['name'],int(n))
	res_warm,res_cool,res_best,res_worst,res_random = [],[],[],[],[]
	for sp in range(int(n)):
		ext_3,new_g=ext_casc(nets_1[p],ext_seq_warm[:sp],preds_0[p],preys_0[p],spp_dict,tre=tre_val,rew=rew_val)	
		res_warm.append(len(ext_3)/n)
	for sp in range(int(n)):
		ext_3,new_g=ext_casc(nets_1[p],ext_seq_cool[:sp],preds_0[p],preys_0[p],spp_dict,tre=tre_val,rew=rew_val)			
		res_cool.append(len(ext_3)/n)
	for sp in range(int(n)):
		ext_3,new_g=ext_casc(nets_1[p],ext_seq_best[:sp],preds_0[p],preys_0[p],spp_dict,tre=tre_val,rew=rew_val)			
		res_best.append(len(ext_3)/n)
	for sp in range(int(n)):
		ext_3,new_g=ext_casc(nets_1[p],ext_seq_worst[:sp],preds_0[p],preys_0[p],spp_dict,tre=tre_val,rew=rew_val)			
		res_worst.append(len(ext_3)/n)
	for sp in range(int(n)):
		ext_3,new_g=ext_casc(nets_1[p],ext_seq_random[:sp],preds_0[p],preys_0[p],spp_dict,tre=tre_val,rew=rew_val)			
		res_random.append(len(ext_3)/n)
	for row in range(int(n)):
		out.write(','.join(map(str,[row/n,res_warm[row],res_cool[row],res_best[row],res_worst[row],res_random[row]]))+'\n')
	print p


out.close()




























