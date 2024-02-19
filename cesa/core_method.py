import numpy as np
from scipy.stats import binom


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
    IP=np.dot(yi,yj)
    #Synchronization
    S_pval=binom.sf((IP+L)/2-1,L,1/2)
    #Anti-synchronization
    AS_pval=binom.cdf((IP+L)/2,L,1/2)
    return S_pval,AS_pval


def make_correction_mat(data,total_read):
    """
    Quantify the dependence of each pair when correcting for relative amounts and create a matrix. 
    The return matrix represents the magnitude of the effect on the species in column j to row i.

    Args:
        data (array_like): Time series matrix with species in rows and time in columns
        total_read (int): Total read count

    Returns:
        correct_mat (np.ndarray): Number of changes in the increase or decrease of species in row i when the variation of species in column j is hypothetically eliminated
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
