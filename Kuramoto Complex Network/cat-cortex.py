import numpy as np
import matplotlib.pyplot as plt
import networkx as nx 
import matplotlib.cm as cm

cat_cortex = np.loadtxt("cat-cortex.csv", unpack="True")
G = nx.DiGraph(cat_cortex)
A = nx.to_dict_of_lists(G,G.nodes())

M = cat_cortex

node_size_in = []
node_size_out = []
color_map = []

cat_cortex_id = ['17','18','19','PLLS','PMLS','AMLS','ALLS','VLS','DLS','21a','21b','20a','20b','7','AES','PS','AI','AII','AAF','P','VP','EPp',
'Tem','3a','3b','1','2','SII','SIV','4g','4','6I','6m','5Am','5AI','5Bm','5BI','SSSAi','SSAo','PFCMil','PFCMd','PFCL','Ia','Ig','CGa','CGp','RS',
'35','36','pSb','Sb','Enr','Hipp']

pos = {0: np.array([ 0.71924299, -0.52181599]), 1: np.array([ 0.54321635, -0.59859512]), 2: np.array([ 0.3483565 , -0.40329295]), 3: np.array([ 0.35666898, -0.17353535]), 4: np.array([ 0.52542354, -0.3143125 ]), 5: np.array([ 0.09517063, -0.46254144]), 6: np.array([ 0.29396561, -0.24544065]), 7: np.array([ 0.85336233, -0.21803655]), 8: np.array([0.56058834, 0.05275141]), 9: np.array([ 0.25522028, -0.45569016]), 10: np.array([ 0.44515112, -0.35899088]), 11: np.array([ 0.45337608, -0.11349625]), 12: np.array([0.24600982, 0.04806786]), 13: np.array([ 0.05634807, -0.20202654]), 14: np.array([-0.01375985, -0.20627626]), 15: np.array([ 0.31368 , -0.058603]), 16: np.array([-0.00402008,  1.        ]), 17: np.array([-0.04141039,  0.61008996]), 18: np.array([-0.06723241,  0.79982502]), 19: np.array([0.24658595, 0.66794997]), 20: np.array([0.1439076 , 0.86297939]), 21: np.array([0.20087867, 0.25891879]), 22: np.array([-0.23343166,  0.55515662]), 23: np.array([-0.20693 , 0.844993]), 24: np.array([-0.61751555, -0.52489754]), 25: np.array([-0.72447873, -0.42511453]), 26: np.array([-0.50288988, -0.47125504]), 27: np.array([-0.46434188, -0.33383946]), 28: np.array([-0.41631963, -0.05817278]), 29: np.array([-0.55626887, -0.32697052]), 30: np.array([-0.47973492, -0.17997691]), 31: np.array([-0.21163196, -0.09363569]), 32: np.array([-0.25355795,  0.04543788]), 33: np.array([-0.59456936, -0.15830478]), 34: np.array([-0.18755204, -0.24683747]), 35: np.array([-0.28585715, -0.3892418 ]), 36: np.array([-0.14848357, -0.40222607]), 37: np.array([-0.34567006, -0.52505567]), 38: np.array([-0.36599892, -0.30572859]), 39: np.array([-0.34207961,  0.46497818]), 40: np.array([-0.16155336,  0.41180297]), 41: np.array([0.02394941, 0.29310196]), 42: np.array([-0.16105092,  0.16196049]), 43: np.array([-0.0849885 ,  0.00625134]), 44: np.array([-0.35457629,  0.12983366]), 45: np.array([0.02554568, 0.11675587]), 46: np.array([0.29015543, 0.3604571 ]), 47: np.array([-0.01765899,  0.19013625]), 48: np.array([-0.02722746,  0.06906645]), 49: np.array([0.38204391, 0.42713071]), 50: np.array([0.06242635, 0.51191121]), 51: np.array([0.21296172, 0.48687901]), 52: np.array([0.53635503, 0.83853754])}


"""
This portion of the code generates the adjacency matrix without weights. 
for i in range(len(cat_cortex)):
    for j in range(len(cat_cortex)):
      if (M[i,j] != 0):
          M[i,j] = 1

np.savetxt("cat-adj.csv",M,fmt="%d")

for i in range(0,19):
    color_map.append('skyblue')
for i in range(19,29):
    color_map.append('red')
for i in range(29,48):
    color_map.append('limegreen')
for i in range(48,65):
    color_map.append('magenta')
"""

for i in range(0,17):
    color_map.append('skyblue')
for i in range(17,24):
    color_map.append('red')
for i in range(24,40):
    color_map.append('limegreen')
for i in range(40,53):
    color_map.append('magenta')

degrees = [val*10 for (node, val) in G.degree()]

ent = dict(G.in_degree(G.nodes()))

sai = dict(G.out_degree(G.nodes()))

for i in range(len(degrees)):
    node_size_out.append(sai[i]*1)

for i in range(len(degrees)):
    node_size_in.append(ent[i]*1)

average_in = sum(node_size_in)/len(node_size_in)-0.5
average_out = sum(node_size_out)/len(node_size_out)

t = np.arange(0,53,1)
average = np.ones_like(t)*average_in

#pos = nx.spring_layout(G)

#setting bars' width
barWidth = 0.325

plt.figure(figsize=(10,5))

p1 = np.arange(len(ent))
p2 = [x + barWidth for x in p1]

plt.bar(p1, node_size_in, color='green', width=barWidth, label= "In degree")
plt.bar(p2, node_size_out, color='purple', width=barWidth, label= "Out degree")

#plt.plot(t,average,color='red',label="average degree")

plt.xlabel("Nodes")
plt.xticks(rotation=90)
plt.ylabel("Degree")
plt.xlim(-0.5,53)
plt.legend()
plt.show()

cat_cortex_id = {0: '$17$',1:'18',2:'19',3:'PLLS',4:'PMLS',5:'AMLS',6:'ALLS',7:'VLS',8:'DLS',9:'21a',10:'21b',11:'20a',12:'20b',13:'7',14:'AES',15:'PS',16:'AI',17:'AII',18:'AAF',19:'P',20:'VP',21:'EPp',22:'Tem',23:'3a',24:'3b',25:'1',26:'2',27:'SII',28:'SIV',29:'4g',30:'4',31:'6I',32:'6m',33:'5Am',34:'5AI',35:'5Bm',36:'5BI',37:'SSSAi',38:'SSAo',39:'PFCMil',40:'PFCMd',41:'PFCL',42:'Ia',43:'Ig',44:'CGa',45:'CGp',46:'RS',
47:'35',48:'36',49:'pSb',50:'Sb',51:'Enr',52:'Hipp'}

nx.draw(G,pos=pos, node_size= degrees, node_color = color_map, with_labels=cat_cortex_id, width=0.15 ,alpha=0.75)
plt.show()


#nx.draw(G,pos=pos, node_size=node_size_out,node_color = color_map, with_labels=True,width=0.15,alpha=0.75)
#plt.show()

#nx.draw(G,pos=pos,node_size=node_size_in,node_color = color_map, labels=ent,width=0.15,alpha=0.75)
#plt.show()


