import sys  
import Data_excluded
import numpy as np
import pandas as pd
import glob2
import glob, os
from io import BytesIO
import datetime as dt
import calendar
from pandas.errors import EmptyDataError
from scipy.ndimage.interpolation import shift


# Function for keeping data that should be used (inbetween starttime and endtime, fog events)
def keepData(data=None,starttime_exclude=None,endtime_exclude=None):
    for start in starttime_exclude.index:
        data_n = data.loc[(starttime_exclude[start]<=data.index)&(data.index<=endtime_exclude[start])]
        if start==starttime_exclude.index[0]:
            data_f=data_n
        else:
            data_f=pd.concat([data_f,data_n],axis=0) #data_f.append(data_n)
    return data_f


#########################################################
### calculating log-norm functions with 1,2,3,4 modes ###
#########################################################
def lognorm(x, N, sigma, xmod):
    return N / (np.sqrt(math.pi * 2) * np.log10(sigma)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod),2) / np.power(np.log10(sigma),2))  

def lognorm2(x, N_s, sigma_s, xmod_s, N, sigma, xmod):
    return (N_s / (np.sqrt(math.pi * 2) * np.log10(sigma_s)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_s),2) / np.power(np.log10(sigma_s),2))) + (N / (np.sqrt(math.pi * 2) * np.log10(sigma)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod),2) / np.power(np.log10(sigma),2))) 

def lognorm3(x, N_vs, sigma_vs, xmod_vs, N_s, sigma_s, xmod_s, N, sigma, xmod):
    return (N_vs / (np.sqrt(math.pi * 2) * np.log10(sigma_vs)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_vs),2) / np.power(np.log10(sigma_vs),2))) + (N_s / (np.sqrt(math.pi * 2) * np.log10(sigma_s)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_s),2) / np.power(np.log10(sigma_s),2))) + (N / (np.sqrt(math.pi * 2) * np.log10(sigma)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod),2) / np.power(np.log10(sigma),2))) 

def lognorm4(x, N_vs, sigma_vs, xmod_vs, N_s, sigma_s, xmod_s, N, sigma, xmod, N_l, sigma_l, xmod_l):
    return (N_vs / (np.sqrt(math.pi * 2) * np.log10(sigma_vs)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_vs),2) / np.power(np.log10(sigma_vs),2))) + (N_s / (np.sqrt(math.pi * 2) * np.log10(sigma_s)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_s),2) / np.power(np.log10(sigma_s),2))) + (N / (np.sqrt(math.pi * 2) * np.log10(sigma)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod),2) / np.power(np.log10(sigma),2))) + (N_l / (np.sqrt(math.pi * 2) * np.log10(sigma_l)) * np.exp( -0.5 * np.power(np.log10(x) - np.log10(xmod_l),2) / np.power(np.log10(sigma_l),2))) 


###############################
### Calculate zero crossing ###
###############################
def zero_crossing(data):
    return np.where(np.diff(np.sign(np.array(data))))


############################################################################
### Function to find the value in a list (items) nearest a value (pivot) ###
############################################################################
# OBS! The items list needs to be sorted from smallest to largest, otherwise it will not work
# upper = True, return closest upper value
# lower = True, return closest lower value
# twosided = TRUE, return the value closest to the pivot regardless of lower or upper
def nearest(items, pivot, upper = False, lower = False, twosided = False):
    
    # Make sure that the items list is sorted in ascending order
    tempdict = {} # Dictionary is used to keep track of the original index (which are used as keys)
    
    for i in range(0, len(items)):
        tempdict[i] = items[i]
    
    tempdict_sorted = dict(sorted(tempdict.items(), key = lambda x:x[1]))
    
    items = list(tempdict_sorted.values()) # new sorted items list
    
    # Find the closest value    
    if pivot > items[0]:
        i1 = max(item for item in items if item < pivot)
    else:
        i1 = items[0]
    
    if pivot < items[-1]:
        i2 = min(item for item in items if item > pivot)
    else:
        i2 = items[-1]

    
    # Orignial indices in the tempdict
    i1_tempdict = np.where(list(tempdict_sorted.values()) == i1)[0][0]
    i2_tempdict = np.where(list(tempdict_sorted.values()) == i2)[0][0]
    
    
    # Return statements based on user input 'upper', 'lower' or 'twosided'
    if twosided == True:        
        if abs(pivot - i1) < abs(pivot - i2):
            return list(tempdict_sorted)[i1_tempdict]
        else:
            return list(tempdict_sorted)[i2_tempdict]
    
    elif lower == True:
        return list(tempdict_sorted)[i1_tempdict]
    
    elif upper == True:
        return list(tempdict_sorted)[i2_tempdict]
    
    else:
        print('Need to choose either upper, lower or twosided')

        
###############################################
### calculate dlogDp (given upper bin size) ###
###############################################
def dlogDp(Dmin,D):
    dlogD=np.log10(D)-shift(np.log10(D),1,cval=np.log10(Dmin))
    return dlogD

###########################################
### calculate Dp (given upper bin size) ###
###########################################
def Dp(Dmin,D):
    D=(D+shift(D,1,cval=Dmin))/2
    return D



##########################################
### Calculate size dist for certain rh ###
##########################################
#######################################
### Calculate size dist for rh=100% ###
#######################################
def rh100(Ddry,kappa,T,rh):
    
    nanmask=~np.isnan(kappa)&~np.isnan(T)&~np.isnan(rh)
    kappa=pd.DataFrame(kappa[nanmask]).dropna()
    T=pd.DataFrame(T[nanmask]).dropna()
    rh=pd.DataFrame(rh[nanmask]).dropna()

    diam=pd.DataFrame(np.linspace(10**(-8),50*10**(-6),100000)).transpose()

    Mw= 0.01801428; # molecular weight of water (kg/mol)
    sigmaw=72.8*10**(-3); # surface tension of water at 
    R=8.314472; # gas constant
    rhow=1000; # density of pure water

    T_K=273.15+T
   
    A=4*Mw*sigmaw/(R*T_K*rhow);

    D100=pd.DataFrame(np.nan,index=kappa.index,
                      columns=Ddry.transpose(),dtype='float')
    i=0
    for D in Ddry:
        numerator = diam**3-D**3
#        print(numerator)
        denominator = pd.DataFrame(diam.values**3-(D**3.*(1-kappa)).values,index=kappa.index)
#        print(denominator)
        frac = pd.DataFrame(numerator.values/denominator.values,index=denominator.index)
#        print(frac)
        exp = np.exp(pd.DataFrame(A.values/diam.values,index=A.index))
#        print(exp)    
        SS=((frac*exp)-1+(1.-rh.values/100.))*100;        
#        print(SS)

#         plt.figure()
#         plt.title('D$_{dry}$ = %.1fnm, $\kappa$=%.2f'%(D*10**9,kappa.mean()))
#         plt.plot(diam.transpose(),SS.mean(axis=0))
#         plt.axhline(y=0,ls='--',c='gray')
#         plt.xscale('log')
#         plt.ylim([-0.3,0.3])
#         plt.xlabel('D$_{p,wet}$ (m)');
#         plt.ylabel('Supersaturaion (%)');
        del numerator,denominator,frac,exp
        for j in range(len(SS)):
#            plt.axvline(x=diam[zero_crossing(SS.iloc[j,:])[-1][-2]].values,c='grey')
            if len(zero_crossing(SS.iloc[j,:])[-1])>2:
                D100.iloc[j,i]=diam[zero_crossing(SS.iloc[j,:])[-1][-2]]
            else:
                D100.iloc[j,i]=np.nan
        i=i+1

    return D100


def rh100_1(Ddry,kappa,T,rh):
    
    nanmask=~np.isnan(kappa)&~np.isnan(T)&~np.isnan(rh)
    kappa=pd.DataFrame(kappa[nanmask]).dropna()
    T=pd.DataFrame(T[nanmask]).dropna()
    rh=pd.DataFrame(rh[nanmask]).dropna()

    diam=pd.DataFrame(np.linspace(10**(-8),50*10**(-6),100000)).transpose()

    Mw= 0.01801428; # molecular weight of water (kg/mol)
    sigmaw=72.8*10**(-3); # surface tension of water at 
    R=8.314472; # gas constant
    rhow=1000; # density of pure water

    T_K=273.15+T
   
    A=4*Mw*sigmaw/(R*T_K*rhow);

    D100=pd.DataFrame(np.nan,index=kappa.index,
                      columns=Ddry.transpose(),dtype='float')
    i=0
    for D in Ddry:
        numerator = diam**3-D**3
#        print(numerator)
        denominator = pd.DataFrame(diam.values**3-(D**3.*(1-kappa)).values,index=kappa.index)
#        print(denominator)
        frac = pd.DataFrame(numerator.values/denominator.values,index=denominator.index)
#        print(frac)
        exp = np.exp(pd.DataFrame(A.values/diam.values,index=A.index))
#        print(exp)    
        SS=((frac*exp)-1+(1.-rh.values/100.))*100;        
#        print(SS)

        plt.figure()
        plt.title('D$_{dry}$ = %.1fnm, $\kappa$=%.2f'%(D*10**9,kappa.mean()))
        plt.plot(diam.transpose(),SS.mean(axis=0))
        plt.axhline(y=0,ls='--',c='gray')
        plt.xscale('log')
        plt.ylim([-0.3,0.3])
        plt.xlabel('D$_{p,wet}$ (m)');
        plt.ylabel('Supersaturaion (%)');
        del numerator,denominator,frac,exp
        if (rh.min().min()==rh.max().max()) and (rh.min().min()==100):
            for j in range(len(SS)):
                D100.iloc[j,i]=diam[zero_crossing(SS.iloc[j,:])[-1][-1]]
                plt.axvline(x=diam[zero_crossing(SS.iloc[j,:])[-1][-1]].values,c='grey')                
        elif rh.min().min()>100:
            for j in range(len(SS)):
                if len(zero_crossing(SS.iloc[j,:])[-1])>1:
                    D100.iloc[j,i]=diam[zero_crossing(SS.iloc[j,:])[-1][-2]]
                    plt.axvline(x=diam[zero_crossing(SS.iloc[j,:])[-1][-2]].values,c='grey')
                else:
                    D100.iloc[j,i]=np.nan
        else:
            D100.iloc[j,i]=np.nan
        i=i+1

    return D100


########################
### refractive index ###
########################
import miepython
def n_water(wavelength):
    """
    Refractive index of water at wavelength.

    Equation is from https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon-24.0C

    Args:
        wavelength: wavelength in microns
    Returns:
        index of refraction
    """
    m_squared = 1.0
    m_squared += 5.666959820E-1 / (1.0 - 5.084151894E-3 / wavelength**2)
    m_squared += 1.731900098E-1 / (1.0 - 1.818488474E-2 / wavelength**2)
    m_squared += 2.095951857E-2 / (1.0 - 2.625439472E-2 / wavelength**2)
    m_squared += 1.125228406E-1 / (1.0 - 1.073842352E1 / wavelength**2)
    refractive_index = np.sqrt(m_squared)
    return refractive_index


##################################
### scattering coefficient Mie ###
##################################
def Q_sca(Dp,lambda_range):
    Qsca=pd.DataFrame(np.nan,index=range(len(Dp)),
                       columns=['D','qext','qsca','qback','g'],dtype='float')

    for i in range(len(Dp)):
        diameter = float(Dp[i])*10**6       # microns
        radius = diameter / 2                      # microns
        num = 200                                  # points to plot
#        lambda_range = 0.88 #np.linspace(0.8, 0.96, num)
        ref_index = n_water(lambda_range)
        x = 2 * np.pi * radius / lambda_range
        qext, qsca, qback, g = miepython.mie(ref_index, x)
        Qsca.loc[i,'D']=diameter
        Qsca.loc[i,'qext']=qext
        Qsca.loc[i,'qsca']=qsca
        Qsca.loc[i,'qback']=qback
        Qsca.loc[i,'g']=g
    return Qsca


##################
### visibility ###
##################
def visibility(dlogDp,Dp,dNdlogDp):
    b=(2.0*math.pi*Dp**2/4*dlogDp*(dNdlogDp*10**6)).sum(axis=1)
    visiblty=3.0/b
    return visiblty

def visibility_Mie(dlogDp,Dp,dNdlogDp,lambda_range):
    Qsca=Q_sca(Dp,lambda_range)
    b=(np.array(dlogDp*Qsca['qext'])*math.pi*Dp**2/4*(dNdlogDp*10**6)).sum(axis=1,min_count=1).replace(0, np.nan)
    visiblty=3.0/(b+13.2*10**(-6)) # 3 comes from the visibility sensor calibration; extinction from Rayleigh see Seinfeld and Pandis (at 520nm)
    return visiblty

def visibility_LWC(LWC): # Kunkel 1984
    b=144.7*LWC**0.88
    visiblty=3.0/b*1000 # as b is given in 1/km
    return visiblty


#########################
### stacked bar chart ###
#########################
def stacked_barchart(ax, tick, series, colors):
    neg_sum = series[series < 0].sum()
    pos_sum = 0.
    for i, val in enumerate(series):
        if val < 0:
            neg_sum -= val
            ax.bar(tick, val, bottom = neg_sum, color = colors[i],
                   align = 'center', width = 0.1, 
                   lw = .25, label = series.index[i])
        else:
            ax.bar(tick, val, bottom = pos_sum, align = 'center',
                   width = 0.1, color = colors[i], 
                   label = series.index[i], lw = .25)
            pos_sum += val
    return ax


#####################
### bivariate fit ###
#####################

"""Function for fitting York, 2004, bivariate fit.
Copyright (C) 2019 Mikko Pitkanen
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
def bivariate_fit(xi, yi, dxi, dyi, ri=0.0, b0=1.0, maxIter=1e6):
    """Make a linear bivariate fit to xi, yi data using York et al. (2004).
    This is an implementation of the line fitting algorithm presented in:
    York, D et al., Unified equations for the slope, intercept, and standard
    errors of the best straight line, American Journal of Physics, 2004, 72,
    3, 367-375, doi = 10.1119/1.1632486
    See especially Section III and Table I. The enumerated steps below are
    citations to Section III
    Parameters:
      xi, yi      x and y data points
      dxi, dyi    errors for the data points xi, yi
      ri          correlation coefficient for the weights
      b0          initial guess b
      maxIter     float, maximum allowed number of iterations
    Returns:
      a           y-intercept, y = a + bx
      b           slope
      S           goodness-of-fit estimate
      sigma_a     standard error of a
      sigma_b     standard error of b
    Usage:
    [a, b] = bivariate_fit( xi, yi, dxi, dyi, ri, b0, maxIter)
    """
    # (1) Choose an approximate initial value of b
    b = b0

    # (2) Determine the weights wxi, wyi, for each point.
    wxi = 1.0 / dxi**2.0
    wyi = 1.0 / dyi**2.0

    alphai = (wxi * wyi)**0.5
    b_diff = 999.0

    # tolerance for the fit, when b changes by less than tol for two
    # consecutive iterations, fit is considered found
    tol = 1.0e-8

    # iterate until b changes less than tol
    iIter = 1
    while (abs(b_diff) >= tol) & (iIter <= maxIter):

        b_prev = b

        # (3) Use these weights wxi, wyi to evaluate Wi for each point.
        Wi = (wxi * wyi) / (wxi + b**2.0 * wyi - 2.0*b*ri*alphai)

        # (4) Use the observed points (xi ,yi) and Wi to calculate x_bar and
        # y_bar, from which Ui and Vi , and hence betai can be evaluated for
        # each point
        x_bar = np.sum(Wi * xi) / np.sum(Wi)
        y_bar = np.sum(Wi * yi) / np.sum(Wi)

        Ui = xi - x_bar
        Vi = yi - y_bar

        betai = Wi * (Ui / wyi + b*Vi / wxi - (b*Ui + Vi) * ri / alphai)

        # (5) Use Wi, Ui, Vi, and betai to calculate an improved estimate of b
        b = np.sum(Wi * betai * Vi) / np.sum(Wi * betai * Ui)

        # (6) Use the new b and repeat steps (3), (4), and (5) until successive
        # estimates of b agree within some desired tolerance tol
        b_diff = b - b_prev

        iIter += 1

    # (7) From this final value of b, together with the final x_bar and y_bar,
    # calculate a from
    a = y_bar - b * x_bar

    # Goodness of fit
    S = np.sum(Wi * (yi - b*xi - a)**2.0)

    # (8) For each point (xi, yi), calculate the adjusted values xi_adj
    xi_adj = x_bar + betai

    # (9) Use xi_adj, together with Wi, to calculate xi_adj_bar and thence ui
    xi_adj_bar = np.sum(Wi * xi_adj) / np.sum(Wi)
    ui = xi_adj - xi_adj_bar

    # (10) From Wi , xi_adj_bar and ui, calculate sigma_b, and then sigma_a
    # (the standard uncertainties of the fitted parameters)
    sigma_b = np.sqrt(1.0 / np.sum(Wi * ui**2))
    sigma_a = np.sqrt(1.0 / np.sum(Wi) + xi_adj_bar**2 * sigma_b**2)

    # calculate covariance matrix of b and a (York et al., Section II)
    cov = -xi_adj_bar * sigma_b**2
    # [[var(b), cov], [cov, var(a)]]
    cov_matrix = np.array(
        [[sigma_b**2, cov], [cov, sigma_a**2]])

    if iIter <= maxIter:
        return a, b, S, cov_matrix
    else:
        print("bivariate_fit.py exceeded maximum number of iterations, " +
              "maxIter = {:}".format(maxIter))
        return np.nan, np.nan, np.nan, np.nan
