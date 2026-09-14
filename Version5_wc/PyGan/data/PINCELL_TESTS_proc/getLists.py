# Python3 Script : UTLTools.py
# Author : R. Guasch, adapted from L. Fede's getLists()
# Purpose : handle issues with generating UTL compatible lists for Dragon5 BU calculations
#
#

import numpy as np

def getLists(burnup_points):

################################## test
    if burnup_points=='test':
        ListeBU=[0.0, 15.0, 30.0, 50.0, 75.0, 150.0, 500.0, 1000.0, 2500.0, 5000.0, 7500.0, 10000.0]
        ListeAUTOP=[1000.0, 15000.0]
        ListeCOMPO=[0.0, 15.0, 30.0, 50.0, 75.0, 150.0, 500.0, 1000.0, 2500.0, 5000.0, 7500.0, 10000.0]

################################## UOx
    elif burnup_points=='UOX':
        ListeBU=[0.0, 0.0300, 0.0500, 0.0750, 0.1500, 0.2500, 0.5000, 0.7500, 1.0000, 2.0000, 2.5000,
                3.0000, 3.5000, 4.0000, 4.5000, 5.0000, 5.5000, 6.0000, 6.5000, 7.0000, 7.5000, 8.0000, 8.5000,
                9.0000, 9.5000, 10.0000, 11.0000, 12.0000, 13.0000, 14.0000, 15.0000, 16.0000, 17.0000, 18.0000, 19.0000, 20.0000,
                22.0000, 24.0000, 26.0000, 28.0000, 30.0000, 32.0000, 36.0000, 40.0000, 44.0000, 48.0000, 52.0000, 56.0000, 60.0000]
        ListeBU = convert_Serpent_to_UTL(ListeBU)
        ListeCOMPO = ListeBU
        ListeAUTOP = [5500.0, 17000.0, 34000.0, 61000.0 ]
    elif burnup_points=="Gd":
        ListeBU=[0.0, 0.015, 0.03, 0.05, 0.075, 0.1125, 0.15, 0.2, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.25, 1.5, 1.75,
                2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75,
                8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75, 12, 12.5, 13.0, 13.5, 14.0, 14.5,15.0,
                16.0, 17.0, 18.0, 19.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0, 52.0, 56.0, 60.0]
        ListeBU = convert_Serpent_to_UTL(ListeBU)
        ListeCOMPO = ListeBU
        ListeAUTOP = [5500.0, 17000.0, 34000.0, 61000.0 ]
################################## free
    elif burnup_points=='free':
        ListeBU=[0.0]
        ListeAUTOP=[0.0]
        ListeCOMPO=[0.0]
    print(f"Length of Burnup points list is : {len(ListeBU)}")
    print("burnup_points: ",burnup_points)
    return [ListeBU,ListeAUTOP,ListeCOMPO]


def convert_Serpent_to_UTL(BU_Serpent_list):
    BU_Serpent_list = np.array(BU_Serpent_list)
    BU_Dragon_list = 1000*BU_Serpent_list
    return BU_Dragon_list


