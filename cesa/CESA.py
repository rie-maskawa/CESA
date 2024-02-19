import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import chi2
from . import core_method as cmd

def COEX(data):
    """
    Return p-values for coexistence and exclusion tests, and an index of bias (ad-bc) for the contingency table.
    In the binary time series of OTU-i and OTU-j, there are 4 cases (1,1), (1,0), (0,1) or (0,0) at each observation time step.
    a, b, c, d represent the percentages of (1,1), (1,0), (0,1) and (0,0), respectively.

    Args:
        data (array_like): Time series matrix with species in rows and time in columns

    Returns:
        (np.adarray, np.adarray, np.adarray)
        C_pval_mat (np.adarray): Upper triangular matrix of p-values for the coexisting test
        E_pval_mat (np.adarray): Upper triangular matrix of p-values for the exclusion test
        ADBC_mat (np.adarray): Upper triangular matrix of ad-bc
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


def SYAS(data,total_read,correction=True):
    """
    Return the p-values for synchronization and anti-synchronization tests, and the cross moments(CM) for each pair

    Args:
        data (array_like): Time series matrix with species in rows and time in columns
        total_read (int): Total read count
        correction (bool): Correction(True) or not(False)

    Returns:
        (np.adarray, np.adarray, np.adarray)
        S_pval_mat (np.adarray): Upper triangular matrix of p-values for the synchronization test
        AS_pval_mat (np.adarray): Upper triangular matrix of p-values for the anti-synchronization test
        CM_mat (np.adarray): Upper triangular matrix of CM
    """

    X=np.array(data)
    N=X.shape[0]

    if correction==True:
        # Number of changes in the increase or decrease of species in row i when the variation of species in column j is hypothetically eliminated
        correct_mat=cmd.make_correction_mat(data,total_read)

        SY_pval_mat=np.zeros((N,N))
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

                SY_pval_list=[]
                AS_pval_list=[]
                
                xi=X[i,:]
                xj=X[j,:]
                #Excluding periods of zero continuity
                indi=(np.diff(xi)!=0) | (xi[:-1]!=0)
                indj=(np.diff(xj)!=0) | (xj[:-1]!=0)
                ind=indi & indj
                L=sum(ind)
                yi=cmd.ternary_transform(xi)[ind]
                yj=cmd.ternary_transform(xj)[ind]
                if L==0:
                    CM_mat[i,j]=0
                else:
                    CM_mat[i,j]=np.dot(yi,yj)/L
                
                #no correction
                SY_pval,AS_pval=cmd.calc_pval(xi,xj)
                SY_pval_list.append(SY_pval)
                AS_pval_list.append(AS_pval)

                #two-body correction
                if ((correct_mat[i,j]!=0) | (correct_mat[j,i]!=0)):
                    #correction
                    xi0=total_read/(total_read-X[j,:])*X[i,:]
                    xj0=total_read/(total_read-X[i,:])*X[j,:]
                    xi=np.where(xi0-X[i,:]>=1,np.floor(xi0),X[i,:])
                    xj=np.where(xj0-X[j,:]>=1,np.floor(xj0),X[j,:])
                    SY_pval,AS_pval=cmd.calc_pval(xi,xj)
                    SY_pval_list.append(SY_pval)
                    AS_pval_list.append(AS_pval)
                #three-body correction
                if len(klist2)!=0:
                    for k in klist2:
                        #correction
                        xi0=total_read/(total_read-X[k,:])*X[i,:]
                        xj0=total_read/(total_read-X[k,:])*X[j,:]
                        xi=np.where(xi0-X[i,:]>=1,np.floor(xi0),X[i,:])
                        xj=np.where(xj0-X[j,:]>=1,np.floor(xj0),X[j,:])
                        SY_pval,AS_pval=cmd.calc_pval(xi,xj)
                        SY_pval_list.append(SY_pval)
                        AS_pval_list.append(AS_pval)

                #Select maximum p-value
                SY_pval_mat[i,j]=np.max(SY_pval_list)
                AS_pval_mat[i,j]=np.max(AS_pval_list)
    
    #no correction
    else:
        SY_pval_mat=np.zeros((N,N))
        AS_pval_mat=np.zeros((N,N))
        CM_mat=np.zeros((N,N))
        for i in range(N):
            for j in range(i+1,N):
                xi=X[i,:]
                xj=X[j,:]
                #Excluding periods of zero continuity
                indi=(np.diff(xi)!=0) | (xi[:-1]!=0)
                indj=(np.diff(xj)!=0) | (xj[:-1]!=0)
                ind=indi & indj
                L=sum(ind)
                yi=cmd.ternary_transform(xi)[ind]
                yj=cmd.ternary_transform(xj)[ind]
                if L==0:
                    CM_mat[i,j]=0
                else:
                    CM_mat[i,j]=np.dot(yi,yj)/L
                
                SY_pval,AS_pval=cmd.calc_pval(xi,xj)

                SY_pval_mat[i,j]=SY_pval
                AS_pval_mat[i,j]=AS_pval

    return SY_pval_mat,AS_pval_mat,CM_mat


def combining_pval(pval_list):
    """
    Returns the matrix of combined p-values for the S samples

    Args:
        pval_list (list): List of p-value matrices

    Returns:
        com_pval (np.ndarray): the matrix of combined p-values for the S samples
    """
    S=len(pval_list)
    Y=-2*np.log(pval_list).sum(axis=0)
    com_pval=chi2.sf(Y,2*S)
    return com_pval
