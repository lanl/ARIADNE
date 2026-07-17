# This module contains pre-defined shapes for correlation matrices

from numpy import diag, ones, zeros, exp, isscalar, shape
from matplotlib import pyplot as plt

def identify_corshape(dim_cor,cor_type,cor_type_arg = {},no_unc = 0):
    """
    Calls functions to generate covariance according to
    the cor_type.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix (no. 
    of rows and columns)
    cor_type : type of correlation matrix
    cor_type_arg : arguments to give for the 
    corrlation shape function
    no_unc : number of uncertainty term, needed 
    to identify cor_type_arg
    ---------------------------------------------

    Output parametersL
    ---------------------------------------------
    cor : correlation matrix
    ---------------------------------------------
    """

    if (dim_cor == 1):
        cor = 1
    else:
        if cor_type == 'Diagonal':
            cor = cor_diagonal(dim_cor)
        elif cor_type == 'Positive_fully':
            cor = cor_positive_fully(dim_cor)
        elif cor_type == 'Constant':
            cor = cor_constant(dim_cor,cor_type_arg,no_unc)
        elif cor_type == 'Gaussian':
            cor = cor_gaussian(dim_cor,cor_type_arg,no_unc)
        elif cor_type == 'Gaussian-Anticorrelated':
            cor = cor_gaussian_anticorrelated(dim_cor,cor_type_arg,no_unc)
        elif cor_type == 'User_defined':
            cor = cor_user_defined(dim_cor,cor_type_arg,no_unc)
        else:
            msg = "Error: Shape of correlation-matrix with type " + cor_type + " not defined."
            raise Exception(msg)              

    return cor

def cor_user_defined(dim_cor,cor_type_arg,no_unc):
    """
    Returns a user-defined correlation matrix.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix,
        will be counter-checked
    cor_type_arg : must contain dictionary 
        "cor_matrix" with number of unc. component
        identifying it.
    no_unc : number of uncertainty component to
    identify the correlation matrix
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    cor : correlation matrix
    ---------------------------------------------
    
    """
    # reading and counterchecking dimension -----
    cor = cor_type_arg['cor_matrix'][str(no_unc)]
    
    if shape(cor)[0] != dim_cor or shape(cor)[1] != dim_cor:
        msg = 'Error: Expecting square matrix of dimension ' +  str(dim_cor) +'x' + str(dim_cor)+'. Got instead'+\
        str(shape(cor)) + ' dimension.'
        raise Exception(msg)
    # --------------------------------------------
    
    # plotting ----------------------------------
    plt.minorticks_on()
    plt.gca().set_aspect('equal')
    plt.title('User-defined correlation for unc. source #'+str(no_unc))
    plt.xlabel('# of data points')
    plt.ylabel('# of data points')
    plt.pcolormesh(cor,cmap='seismic',vmin=-1,vmax=1)
    plt.colorbar()
    plt.show()
    # --------------------------------------------
    
    return cor

def cor_gaussian_anticorrelated(dim_cor,cor_type_arg,no_unc):
    """
    Gives anticorrelated matrix with a certain
    turning-point energy softened by a Gaussian.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix
    cor_type_arg : must contain outgoing neutron
          energy and turning-point energy of
          the anti-correlated matrix, and 
          may contain an additional parameter 
          "damp_term < 1"
    no_unc : number of uncertainty component to
    identify the damp_term
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    cor : correlation matrix
    ---------------------------------------------
    
    """
    cor = zeros([dim_cor,dim_cor],dtype = float)

    # anti-correlated  matrix -------------------    
    if 'eout' in cor_type_arg:
        E = cor_type_arg['eout']
    if 'einc' in cor_type_arg:
        E = cor_type_arg['einc']

    if isscalar(cor_type_arg['eout-turningpoint']):
        Eturn = cor_type_arg['eout-turningpoint']
    else:
        Eturn = cor_type_arg['eout-turningpoint'][no_unc]

        

    for index1 in range(0,dim_cor):
        for index2 in range(0,dim_cor):
            if (E[index1] > Eturn) and (E[index2] < Eturn) :
                cor[index1,index2] = -1
            elif (E[index1] < Eturn) and (E[index2] > Eturn) :
                cor[index1,index2] = -1
            else:
                cor[index1,index2] = 1                
    # -------------------------------------------

    # Gaussian softening ------------------------
    if 'damp_term' in cor_type_arg:
        if isscalar(cor_type_arg['damp_term']):
            damp_term = cor_type_arg['damp_term']
        else:
            damp_term = cor_type_arg['damp_term'][no_unc]
    else:
        damp_term = 1.0

    if (damp_term > 1.0 ) or (damp_term < 0.0):
        msg = "Warning: Damp_term of Gaussian correlation matrix shape can only assume value [0,1] but values is" + str(damp_term)+", value set to 1.0"
        print(msg)
        damp_term = 1.0

    for index1 in range(0,dim_cor):
        for index2 in range(0,dim_cor):
            dscr = ((E[index1]-E[index2])*damp_term/(max(E[index1],E[index2])))**2.0
            cor[index1,index2] = cor[index1,index2]*exp(-1.0*dscr)
    # -------------------------------------------

    return cor


def cor_diagonal(dim_cor):
    """
    Gives diagonal correlation matrix.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    diagonal correlation matrix
    ---------------------------------------------
    """
    return diag(ones(dim_cor,dtype = float))



def cor_gaussian(dim_cor,cor_type_arg,no_unc):
    """
    Gives gaussian-curve shaped correlation matrix.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix
    cor_type_arg : must contain outgoing neutron
          energy, may contain an additional 
          parameter "damp_term < 1"
    no_unc : number of uncertainty component to
    identify the damp_term
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    cor : correlation matrix
    ---------------------------------------------
    """
    cor = zeros([dim_cor,dim_cor],dtype = float)
    
    if 'eout' in cor_type_arg:
        E = cor_type_arg['eout']
    if 'einc' in cor_type_arg:
        E = cor_type_arg['einc']

    if 'damp_term' in cor_type_arg:
        if isscalar(cor_type_arg['damp_term']):
            damp_term = cor_type_arg['damp_term']
        else:
            damp_term = cor_type_arg['damp_term'][no_unc]
    else:
        damp_term = 1.0

    if (damp_term > 1.0 ) or (damp_term < 0.0):
        msg = "Warning: Damp_term of Gaussian correlation matrix shape can only assume value [0,1] but values is" + str(damp_term)+", value set to 1.0"
        print(msg)
        damp_term = 1.0

    for index1 in range(0,dim_cor):
        for index2 in range(0,dim_cor):
            dscr = ((E[index1]-E[index2])*damp_term/(max(E[index1],E[index2])))**2.0
            cor[index1,index2] = exp(-1.0*dscr)

    return cor

def cor_positive_fully(dim_cor):
    """
    Gives fully positively correlated correlation 
    matrix.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    fully positively correlated correlation matrix
    ---------------------------------------------
    """
    return ones([dim_cor,dim_cor],dtype = float)


def cor_constant(dim_cor,cor_type_arg,no_unc):
    """
    Gives positively correlated correlation 
    matrix which is constant in its off-diagonal
    entries.

    Input parameters:
    ---------------------------------------------
    dim_cor : dimension of correlation matrix
    
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    fully positively correlated correlation matrix
    cor_type_arg : must contain an parameter 
    "damp_term <= 1" and "damp_term >= 0" which
    corresponds to the off-diagonal, constant,
    correlation factor.
    no_unc : number of uncertainty component to
    identify the damp_term
    ---------------------------------------------
    """

    # assignment of off-diagonal correlation factor (damp_term) 
    if 'damp_term' in cor_type_arg:
        if isscalar(cor_type_arg['damp_term']):
            damp_term = cor_type_arg['damp_term']
        else:
            damp_term = cor_type_arg['damp_term'][no_unc]
    else:
        damp_term = 1.0

    if (damp_term > 1.0 ) or (damp_term < 0.0):
        msg = "Warning: Damp_term of Constant correlation matrix shape can only assume value [0,1] but values is" + str(damp_term)+", value set to 1.0"
        print(msg)
        damp_term = 1.0
    # ----------------------------------------------

    # assignment of correlation matrix -------------
    cor = ones([dim_cor,dim_cor],dtype = float)*damp_term 
    for index1 in range(0,dim_cor):
        cor[index1,index1] = 1.0
    # ----------------------------------------------
    
    return cor
