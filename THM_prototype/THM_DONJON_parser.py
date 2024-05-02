# THM_DONJON_parser class
# Used to post treat THM_DONJON data and compare it to THM_prototype results
# Authors : Clement Huet, Raphael Guasch
import numpy as np
class THM_DONJON_parser:
    def __init__(self, path, t0, time_step, t_end, nz, h_mesh):
        self.path = path
        self.data = None
        self.lines = self.load_data()
        
        # Axial mesh information
        self.nz = nz
        print(f"self.nz is {self.nz}")
        self.h_mesh = h_mesh
        self.dz = self.h_mesh/self.nz
        

        # Time mesh information
        self.begin_time = t0
        self.dt = time_step
        self.end_time = t_end
        if self.dt == 0:
            self.mumber_time_steps = 1
        else:
            self.mumber_time_steps = int((self.end_time-self.begin_time) / self.dt)
   
        self.list_tables = self.parse_lines()
        self.global_table = self.create_global()
        self.TCOMB, self.TSURF, self.DCOOL, self.TCOOL, self.PCOOL, self.HCOOL, self.QFUEL, self.QCOOL, self.VOID, self.QUAL, self.SLIP, self.FLOW_REGIME = self.create_values()

    def load_data(self):
        with open(self.path, 'r') as data_file:
            lines=data_file.readlines()
        return lines


    def parse_lines(self):
        #Parse the lines from the input file to create a list of arrays
        #Each sub-array contains the time at the 0th index and the values (list) of fields (values) at the 1st index
        t=0
        i=0
        table_time_val=[]
        while i<len(self.lines):
            if " ________________________________________________________________________________________________________________________________________________________________" in self.lines[i]:
                #print(f"Itération trouvée a la ligne i = {i}")
                #print(f'Le pas de temps vaut t = {t} s')
                table_time_val.append([t,self.lines[i+4:i+self.nz+4]])
                t+=1
                i+=self.nz
            else:
                i+=1
        return table_time_val


    def parse_table(self, table_t):
        # takes the index of the list created by parse_lines as input
        # returns the associated time step and processed array
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
    

    def create_global(self):
        # create the global array of results which concatenates all the arrays for
        # each time step : data needed to create the values with the create_values function
        big_table = []
        for table in self.list_tables:
            t,tableau = self.parse_table(table)
            big_table.append(tableau)

        return big_table


    def create_values(self):
        # create arrays with the values of the fields for each time step
        # Lines are the time steps, columns are the vertical mesh elements
        # creates the list of time steps and the list of vertical steps
        if self.dt != 0:
            T=np.arange(0,self.end_time,self.dt)
        else:
            T=np.array([1])
        list_z=np.arange(0,self.h_mesh,self.dz)
        time = len(self.global_table)
        TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME=np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz)),np.zeros((time, self.nz))
        for i in range(time):
            for j in range(self.nz):
                TCOMB[i][j]=self.global_table[i][j][0]
                TSURF[i][j]=self.global_table[i][j][1]
                DCOOL[i][j]=self.global_table[i][j][2]
                TCOOL[i][j]=self.global_table[i][j][3]
                PCOOL[i][j]=self.global_table[i][j][4]
                HCOOL[i][j]=self.global_table[i][j][5]
                QFUEL[i][j]=self.global_table[i][j][6]
                QCOOL[i][j]=self.global_table[i][j][7]
                VOID[i][j]=self.global_table[i][j][8]
                QUAL[i][j]=self.global_table[i][j][9]
                SLIP[i][j]=self.global_table[i][j][10]
                FLOW_REGIME=self.global_table[i][j][11]
        return TCOMB, TSURF, DCOOL, TCOOL, PCOOL, HCOOL, QFUEL, QCOOL, VOID, QUAL, SLIP, FLOW_REGIME





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