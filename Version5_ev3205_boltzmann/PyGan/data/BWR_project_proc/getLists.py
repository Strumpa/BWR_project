# Python3 Srcipt : UTLTools.py
# Author : R. Guasch, adapted from L. Fede's getLists()
# Purpose : handle issues with generating UTL compatible lists for Dragon5 BU calculations
#


#

def getLists(burnup_points):
    
#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Definition of ListBU ListAUTOP ListCOMPO            #
#                                                                   #
#####################################################################

################################## test
    if burnup_points=='test':
        ListeBU=[0.0,    15.00]
        ListeAUTOP=[15.00]
        ListeCOMPO=[0.0,    15.00]

################################## UOx
    elif burnup_points=='UOx':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeAUTOP=[30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='UOx_autop5':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeAUTOP=[5500.0, 17000.0, 34000.0, 61000.0  ]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='UOx2':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0,  28000.0,  30000.0, 32000.0, 34000.0, 36000.0, 38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0, 50000.0, 52000.0, 54000.0, 56000.0, 58000.0, 60000.0]
        ListeAUTOP=[30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0,  28000.0,  30000.0, 32000.0, 34000.0, 36000.0, 38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0, 50000.0, 52000.0, 54000.0, 56000.0, 58000.0, 60000.0]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='UOx2_autop5':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0,  28000.0,  30000.0, 32000.0, 34000.0, 36000.0, 38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0, 50000.0, 52000.0, 54000.0, 56000.0, 58000.0, 60000.0]
        ListeAUTOP=[5500.0, 17000.0, 34000.0, 61000.0  ]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='UOx4':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0, 25000.0, 26000.0, 27000.0,  28000.0, 29000.0, 30000.0, 31000.0, 32000.0, 33000.0, 34000.0, 35000.0, 36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 56000.0, 57000.0, 58000.0, 59000.0, 60000.0]
        ListeAUTOP=[30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0, 25000.0, 26000.0, 27000.0,  28000.0, 29000.0, 30000.0, 31000.0, 32000.0, 33000.0, 34000.0, 35000.0, 36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 56000.0, 57000.0, 58000.0, 59000.0, 60000.0]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='UOx4_autop5':
        ListeBU=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0, 23000.0, 24000.0, 25000.0, 26000.0, 27000.0,  28000.0, 29000.0, 30000.0, 31000.0, 32000.0, 33000.0, 34000.0, 35000.0, 36000.0, 37000.0, 38000.0, 39000.0, 40000.0, 41000.0, 42000.0, 43000.0, 44000.0, 45000.0, 46000.0, 47000.0, 48000.0, 49000.0, 50000.0, 51000.0, 52000.0, 53000.0, 54000.0, 55000.0, 56000.0, 57000.0, 58000.0, 59000.0, 60000.0]
        ListeAUTOP=[5500.0, 17000.0, 34000.0, 61000.0  ]
        ListeCOMPO=[0.0, 30.0, 50.0, 75.0, 150.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]


################################## Gadolinium
    elif burnup_points=='Gd':
        ListeBU=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeAUTOP=[15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeCOMPO=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='Gd_autop3':
        ListeBU=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeAUTOP=[1000.0, 10000.0, 20000.0,  65000.0]
        ListeCOMPO=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
    elif burnup_points=='Gd_autop4':
        ListeBU=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]
        ListeAUTOP=[15000.0, 65000.0]
        ListeCOMPO=[0.0, 15.0, 30.0, 50.0, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0, 625.0, 750.0, 875.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0, 11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 24000.0, 28000.0, 32000.0, 36000.0, 40000.0, 44000.0, 48000.0, 52000.0, 56000.0, 60000.0]

################################## CASMO-5 BU steps for comparison

    elif burnup_points=="CASMO-5":
        ListeBU = [0., 100., 500., 1000., 1500., 2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.,  7500.,  8000.,  8500.,  9000.,  9500., 10000., 10500., 11000., 11500., 12000., 12500., 13000., 13500., 14000., 14500., 15000., 15500., 16000., 16500., 17500., 20000., 22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000.]
        ListeAUTOP = [0., 100., 500.,  1000.,  1500.,  2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.,  7500.,  8000.,  8500.,  9000.,  9500., 10000., 10500., 11000., 11500., 12000., 12500., 13000., 13500., 14000., 14500., 15000., 15500., 16000., 16500., 17500., 20000., 22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000., 41000.]
        ListeCOMPO = [0., 100., 500.,  1000.,  1500.,  2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.,  7500.,  8000.,  8500.,  9000.,  9500., 10000., 10500., 11000., 11500., 12000., 12500., 13000., 13500., 14000., 14500., 15000., 15500., 16000., 16500., 17500., 20000., 22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000.]
    elif burnup_points=="CASMO-5_autop5":
        ListeBU = [0., 100.,   500.,  1000.,  1500.,  2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.,  7500.,  8000.,  8500.,  9000.,  9500., 10000., 10500., 11000., 11500., 12000., 12500., 13000., 13500., 14000., 14500., 15000., 15500., 16000., 16500., 17500., 20000., 22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000.]
        ListeAUTOP = [2000., 5000., 14500., 27500., 35000., 41000. ]
        ListeCOMPO = [0., 100., 500.,  1000.,  1500.,  2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.,  7500.,  8000.,  8500.,  9000.,  9500., 10000., 10500., 11000., 11500., 12000., 12500., 13000., 13500., 14000., 14500., 15000., 15500., 16000., 16500., 17500., 20000., 22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000.]

################################## free
    elif burnup_points=='free':
        ListeBU=[0.0]
        ListeAUTOP=[0.0]
        ListeCOMPO=[0.0]
    print(f"Length of Burnup points list is : {len(ListeBU)}")
    return [ListeBU,ListeAUTOP,ListeCOMPO]


