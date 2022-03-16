import numpy as np
from scipy.optimize import fsolve
TI = 650;
echos = 5;
n = 30;
HR = 65;
TR = 1/HR*60*1000;
flipA = 15;
M0 = 1;
PMR = 1.3477;
T1_myo = 1229;
#T1_pl = 900;
def eq(T1_pl):
    return(abs(M0*(np.exp(-TI/T1_pl) - 1) - (np.exp(-TI/T1_pl)*(M0*(np.exp((echos*n - TR + TI)
/T1_pl) - 1) + np.exp((echos*n - TR + TI)/T1_pl)*(M0*np.exp(-(echos*n)/T1_pl)*np.cos((np.pi*flipA)
/180)**n*(np.exp(-TI/T1_pl) - 1) + (M0*(np.exp(-echos/T1_pl) - 1)*((np.exp(-echos/T1_pl)*np.cos((np.pi*flipA)/180))**n - 1))
/(np.exp(-echos/T1_pl)*np.cos((np.pi*flipA)/180) - 1))))/(np.exp(-(echos*n)/T1_pl)*np.exp(-TI/T1_pl)*np.exp((echos*n - TR + TI)
/T1_pl)*np.cos((np.pi*flipA)/180)**n + 1))/abs(M0*(np.exp(-TI/T1_myo) - 1) - (np.exp(-TI/T1_myo)*(M0*(np.exp((echos*n - TR + TI)
/T1_myo) - 1) + np.exp((echos*n - TR + TI)/T1_myo)*(M0*np.exp(-(echos*n)/T1_myo)*np.cos((np.pi*flipA)/180)**n*(np.exp(-TI/T1_myo)
- 1) + (M0*(np.exp(-echos/T1_myo) - 1)*((np.exp(-echos/T1_myo)*np.cos((np.pi*flipA)/180))**n - 1))/(np.exp(-echos/T1_myo)
*np.cos((np.pi*flipA)/180) - 1))))/(np.exp(-(echos*n)/T1_myo)*np.exp(-TI/T1_myo)*np.exp((echos*n - TR + TI)/T1_myo)
*np.cos((np.pi*flipA)/180)**n + 1))-PMR);

T1_pl =fsolve(eq,1000)

print(T1_pl)