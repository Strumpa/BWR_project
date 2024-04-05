import numpy as np
import matplotlib.pyplot as plt

###########
dt = 0.5
end_time = 20.0
n_z = 20
h_mesh = 80
###########

n_pas_t = int(end_time / dt)
dz = h_mesh/n_z
###########

def lire(nom_fichier):
    with open(nom_fichier, 'r') as fichier:
      lignes=fichier.readlines()
    return lignes

#prend les lignes du fichier en entrée et crée une liste des tableaux
#chaque sous liste à le temps en indice 0 et le tableau liée en indice 2
def sparse_lignes(lignes):
    t=0
    i=0
    table_time_val=[]
    while i<len(lignes):
        if " ________________________________________________________________________________________________________________________________________________________________" in lignes[i]:
            #print(f"Itération trouvée a la ligne i = {i}")
            #print(f'Le pas de temps vaut t = {t} s')
            table_time_val.append([t,lignes[i+4:i+n_z+4]])
            t+=1
            i+=n_z
        else:
            i+=1
    return table_time_val

#prend un indice de la liste crée par sparse_ligne en entrée
#et return le pas de temps et le tableau traité associé
def sparse_table(table_t):
    t = table_t[0]
    table_temp = table_t[1]
    tableau=[]
    for i in range(len(table_temp)):
        table_temp[i] = table_temp[i].strip()
        table_temp[i] = table_temp[i].replace(" ", "")
        table_temp[i] = table_temp[i].replace("|", ",")
        table_temp[i] = table_temp[i].split(',')
        table_temp[i] = table_temp[i][2:-1]
        for j in range(len(table_temp[i])):
            table_temp[i][j]=float(table_temp[i][j])
        tableau.append(table_temp[i])
    return t, tableau

#créer les tableau de chaque valeurs les lignes sont les pas
#de temps et chaque colonne un element de maillage vertical
#créer également la liste des temps et la liste des pas verticaux
def create_values(big_table):
    T=np.arange(0,end_time,dt)
    print(T)
    list_z=np.arange(0,h_mesh,dz)
    time = len(big_table)
    z_mesh = len(big_table[0])
    TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME=np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh)),np.zeros((time, z_mesh))
    for i in range(time):
        for j in range(z_mesh):
            TCOMB[i][j]=big_table[i][j][0]
            TSURF[i][j]=big_table[i][j][1]
            DCOOL[i][j]=big_table[i][j][2]
            TCOOL[i][j]=big_table[i][j][3]
            PCOOL[i][j]=big_table[i][j][4]
            HCOOL[i][j]=big_table[i][j][5]
            QFUEL[i][j]=big_table[i][j][6]
            QCOOL[i][j]=big_table[i][j][7]
            VOID[i][j]=big_table[i][j][8]
            QUAL[i][j]=big_table[i][j][9]
            SLIP[i][j]=big_table[i][j][10]
            FLOW_REGIME=big_table[i][j][11]
    return TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME

#crée le tableau global qui concatene tout les tableau pour 
#chaque valeur de t dont a besoin create_value
def create_global(list_table):
    big_table = []
    for table in list_table:
        t,tableau = sparse_table(table)
        big_table.append(tableau)
    return big_table

#commandes pour recuperer  TCOMB à t=0 puis à x=0
lignes = lire('thm.result')
list_table=sparse_lignes(lignes)
big_table=create_global(list_table)
TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME = create_values(big_table)
#print('TCOMB',TCOMB)

def recup_val_tcst(val, t):
    t=float(t)
    T=list(np.arange(0,end_time,dt))
    list_z=list(np.arange(0,h_mesh,dz))
    return list_z, val[T.index(t)]

def recup_val_xcst(val, x):
    x=float(x)
    list_z=list(np.arange(0,h_mesh,dz))
    T=list(np.arange(0,end_time,dt))
    return T,val[:,list_z.index(x)]

lignes = lire('thm.result')
list_table=sparse_lignes(lignes)
big_table=create_global(list_table)
TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME = create_values(big_table)
#print('TCOMB',TCOMB)
X,TCOMB_Tcst=recup_val_tcst(TCOMB,0)
T,TCOMB_Xcst=recup_val_xcst(TCOMB,0)
print(recup_val_tcst(TCOMB,0))
print(recup_val_xcst(TCOMB,0))

fig, ax = plt.subplots( nrows=1, ncols=1 )
plt.plot(X,TCOMB_Tcst, label="Température combustible à t=0")
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.legend()
fig.savefig('result.png')