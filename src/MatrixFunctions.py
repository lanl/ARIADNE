# This module contains matrix functions

from numpy import linalg as LA
from numpy import allclose, shape

def get_cor(cov):
    """
    Calculates correlation matrix from a covariance matrix.
    
    Input parameters:
    ------------------------------------------------
    cov : two-dimensional array to be converted to
    cor.
    ------------------------------------------------

    Output parameters:
    ------------------------------------------------
    cor : two-dimensional array, correlation matrix.
    ------------------------------------------------
    """
    from numpy import shape, zeros, sqrt
    
    # check matrix ---------------------------------
    # check if two-dimensional .....................
    if len(shape(cov)) != 2:
        msg = 'Error: Expecting 2-dim array, got instead ' + \
        str(len(shape(cov))) + ' dimensions.'
        raise Exception(msg)
    # ..............................................

    # check if square matrix .......................
    if (shape(cov)[0]) != (shape(cov)[1]):
        msg = 'Error: Cov is not square.'
        raise Exception(msg)
    # ..............................................
    # ----------------------------------------------

    # getting cor ----------------------------------
    dim = shape(cov)[0]
    cor = zeros([dim,dim],dtype = float)

    for index1 in range(0,dim):
        for index2 in range(0,dim):        
            cor[index1,index2] = cov[index1,index2]/\
                sqrt(cov[index1,index1]*cov[index2,index2])
    # ----------------------------------------------

    return cor

def checking_cov(cov):
    """
    Tests covariance for a) symmetry and b) positive semi-definitness

    Input parameters:
    ---------------------------------------------
    cov : two-dimensional array of cov to be tested
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    none, error messages are given by the sub-program
    ---------------------------------------------
    """

    if shape(shape(cov))[0] != 0:

        msg = check_symmetry(cov)
        print(msg)

        msg = check_positivesemidefiniteness(cov)
        print(msg)
    else:
        print('Covariance is 1-dimensional. No symmetry and postive semi-definiteness checked.')

    return None


def check_positivesemidefiniteness(cov):
    """
    Tests covariance for positive semi-definitness.

    Input parameters:
    ---------------------------------------------
    cov : two-dimensional array of cov to be tested
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    msg : success, warning or error message
    ---------------------------------------------
    """
    bound = -1.0e10

    eigval = LA.eigvalsh(cov)

    if (all(eigval) >= 0):
        msg = "Ok: Matrix is positive semi-definite."
    else:
        if (all(eigval) >= bound):
            msg = "Warning: Marix has small negative eigenvalues, i.e. smaller than" + str(bound)
        else:
            msg = "Error: Matrix is not positive semi-definite. Smallest negative eigenvalue:" + str(min(eigval))
            raise Exception(msg)

    return msg


def check_symmetry(cov):
    """
    Tests covariance for symmtery.

    Input parameters:
    ---------------------------------------------
    cov : two-dimensional array of cov to be tested
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    msg : success, warning or error message
    ---------------------------------------------
    """

    if allclose(cov.transpose(), cov):
        msg = "Ok: Matrix is symmetric."
    else:
        msg = "Error: Matrix is not symmetric."
        raise Exception(msg)

    return msg
