import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.stats import fisher_exact
from scipy.stats import chi2


def ternary_transform(x):
    """
    Ternary transformation of time series by increase/decrease

    Args:
        x (array_like): Time series

    Returns:
        y (np.ndarray): Ternary transformed time series
    """

    y=np.where(np.diff(x)>0,1,np.diff(x))
    y=np.where(y<0,-1,y)
    return y
 

def calc_pval(xi,xj):
    """
    Calculate p-values for synchronization and anti-synchronization for pair (i,j)

    Args:
        xi (array_like): Time series of species i
        xj (array_like): Time series of species j

    Returns:
        S_pval, AS_pval (np.float64, np.float64) : p-value of synchronization, p-value of anti-synchronization of pair (i,j)
    """

    #Excluding periods of zero continuity
    indi=(np.diff(xi)!=0) | (xi[:-1]!=0)
    indj=(np.diff(xj)!=0) | (xj[:-1]!=0)
    ind=indi & indj
    L=sum(ind)
    yi=ternary_transform(xi)[ind]
    yj=ternary_transform(xj)[ind]
    C=sum(yi*yj==1)
    #Synchronization
    S_pval=binom.sf(C,L,1/2)
    #Anti-synchronization
    AS_pval=binom.cdf(C,L,1/2)
    return S_pval,AS_pval


def make_correction_mat(data,total_read):
    """
    _summary_

    Args:
        df (array_like): _description_
        total_read (int): シークエンスのリード数

    Returns:
        _type_: _description_
    """

    X=np.array(data)
    N=X.shape[0]
    #Stores the number of signs of increase or decrease in species i changed by species j
    correct_mat=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            xi0=total_read/(total_read-X[j,:])*X[i,:]
            xj0=total_read/(total_read-X[i,:])*X[j,:]
            #If the change in correction is greater than 1, the correction value rounded down to the nearest whole number is adopted
            xi=np.where(xi0-X[i,:]>=1,np.floor(xi0),X[i,:])
            xj=np.where(xj0-X[j,:]>=1,np.floor(xj0),X[j,:])
            yi=ternary_transform(xi)
            yj=ternary_transform(xj)
            correct_mat[i,j]=np.sum(yi!=ternary_transform(X[i,:]))
            correct_mat[j,i]=np.sum(yj!=ternary_transform(X[j,:]))
    return correct_mat


def CE(data):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """

    X=np.array(data)
    N=X.shape[0]
    T=X.shape[1]
    Y=np.where(X!=0,1,0)

    C_pval_mat=np.zeros((N,N))
    E_pval_mat=np.zeros((N,N))
    ADBC_mat=np.zeros((N,N))
    for i in range(N):
        yi=Y[i,:]
        for j in range(i+1,N):
            yj=Y[j,:]
            a=np.sum((yi==1) & (yj==1))
            b=np.sum((yi==1) & (yj==0))
            c=np.sum((yi==0) & (yj==1))
            d=np.sum((yi==0) & (yj==0))
            C_pval_mat[i,j]=fisher_exact([[a,b],[c,d]],'greater')[1]
            E_pval_mat[i,j]=fisher_exact([[a,b],[c,d]],'less')[1]
            ADBC_mat[i,j]=(a*d-b*c)/T**2

    return C_pval_mat,E_pval_mat,ADBC_mat


def SA(df,correct_mat,total_read):
    """
    _summary_

    Args:
        df (_type_): _description_
        correct_mat (_type_): _description_
        total_read (_type_): _description_

    Returns:
        _type_: _description_
    """

    correct_mat=np.array(correct_mat)
    N=df.shape[0]
    X=np.array(df)

    S_pval_mat=np.zeros((N,N))
    AS_pval_mat=np.zeros((N,N))
    CM_mat=np.zeros((N,N))
    for i in range(N):
        for j in range(i+1,N):
            #correction
            #If there is a species that causes a sign change in either i or j, calculate p-values for all the corrections by those species and adopt the largest value.
            klist_i,=np.where(correct_mat[i,:]!=0)
            klist_j,=np.where(correct_mat[j,:]!=0)
            klist1=np.array(list(set(set(klist_i) | set(klist_j))))
            klist2=klist1[(klist1!=i) & (klist1!=j)]

            S_pval_list=[]
            AS_pval_list=[]
            
            #no correction
            xi=X[i,:]
            xj=X[j,:]
            S_pval,AS_pval=calc_pval(xi,xj)
            S_pval_list.append(S_pval)
            AS_pval_list.append(AS_pval)

            #two-body correction
            if (i in klist1) | (j in klist1):
                #correction
                xi0=total_read/(total_read-X[j,:])*X[i,:]
                xj0=total_read/(total_read-X[i,:])*X[j,:]
                xi=np.where(xi0-X[i,:]>=1,np.floor(xi0),X[i,:])
                xj=np.where(xj0-X[j,:]>=1,np.floor(xj0),X[j,:])
                S_pval,AS_pval=calc_pval(xi,xj)
                S_pval_list.append(S_pval)
                AS_pval_list.append(AS_pval)
            #three-body correction
            if len(klist2)!=0:
                for k in klist2:
                    #correction
                    xi0=total_read/(total_read-X[k,:])*X[i,:]
                    xj0=total_read/(total_read-X[k,:])*X[j,:]
                    xi=np.where(xi0-X[i,:]>=1,np.floor(xi0),X[i,:])
                    xj=np.where(xj0-X[j,:]>=1,np.floor(xj0),X[j,:])
                    S_pval,AS_pval=calc_pval(xi,xj)
                    S_pval_list.append(S_pval)
                    AS_pval_list.append(AS_pval)

            S_pval_mat[i,j]=np.max(S_pval_list)
            AS_pval_mat[i,j]=np.max(AS_pval_list)

            #Excluding periods of zero continuity
            indi=(np.diff(xi)!=0) | (xi[:-1]!=0)
            indj=(np.diff(xj)!=0) | (xj[:-1]!=0)
            ind=indi & indj
            L=sum(ind)
            yi=ternary_transform(xi)[ind]
            yj=ternary_transform(xj)[ind]
            CM_mat[i,j]=np.dot(yi,yj)/L

    return S_pval_mat,AS_pval_mat,CM_mat


def combining_pval(pval_list):
    """
    _summary_

    Args:
        pval_list (_type_): _description_

    Returns:
        _type_: _description_
    """
    S=len(pval_list)
    Y=-2*np.log(pval_list).sum(axis=0)
    com_pval=chi2.sf(Y,2*S)
    return com_pval
