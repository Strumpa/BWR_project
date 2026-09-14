## Helpers to assist in articulating multiphysics iterations.
# 
# ########## Begin helper functions declaration ##########
import numpy as np
import matplotlib.pyplot as plt
import os

## Fonction used for the convergence algorithm
def underRelaxation(Field, OldField, underRelaxationFactor):
    return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

def convergence_TH(Field, OldField, tol):
    conv=False
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i]) > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_POW(Field, OldField, tol):
    """
    Tolerance in % of the nodal power, to allow for an easier comparison between different total powers.
    """
    conv=False
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i])/OldField[i] > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_pcm(keffNew, keffOld, tol):
    conv=False
    print(f"Testing convegence on keff with residual = {np.abs(keffNew - keffOld)*1e5} pcm")
    if np.abs(keffNew - keffOld)*1e5 > tol:
        conv = False
    else:
        conv = True
    return conv

def compute_RMS(Field):
    return np.sqrt(np.sum(Field**2)/len(Field))

def compute_residuals(field):
    """
    Compute the residuals of the field
    field : list : list of the fields to compute the residuals in %
    """
    residuals = (field[-1] - field[-2])*100/field[-2]
    print(f"Residuals of the field = {residuals} %")
    return residuals

def quickPlot(x, y, title, xlabel, ylabel, saveName, path, SAVE_DIR):
    fig,ax = plt.subplots(figsize=(10,6))
    #print(f"In quick plot save at {SAVE_DIR}, from path = {path}, with save_name = {saveName}")
    if len(y) == len(x) and isinstance(y[0], np.float32):
        ax.scatter(x, y)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    else:
        for i in range(len(y)):
            if i%10==0:
                data = y[i]
                ax.plot(x, data, '2-',linewidth=1, label=f"iteration {i}")
        ax.plot(x, y[-1], '2-',linewidth=1, label=f"iteration {len(y)}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    plt.close(fig)
    return


######## End helper functions declaration ##########
