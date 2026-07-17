# This module contains basic physics function

from numpy import isscalar, array, shape, zeros, arange

def maxw(T,T_unit,x,x_unit): 
    """
    Calculates a Maxwellian distribution with temperature
    T on lattices x.
    
    Input parameters:
    ---------------------------------------------
    T : temperature of Maxwellian
    x : lattice, either one-dimensional value or
    one-dimensional array
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    Maxwellian distribution.
    ---------------------------------------------
    """
    from numpy import shape, exp, sqrt
    from MathsPhysicsConstants import PI

    # test if right dimension -------------------
    if (len(shape(x)) > 1) :
        msg = 'Error: Expect either a scaler or 1 dimension array. \
        Instead, the lattice shape is ' + str(shape(x))
        raise Exception(msg)
    # -------------------------------------------
    
    # test if same unit -------------------------
    if T_unit != x_unit:
        msg = 'Error: Unit of T and lattice of Maxwellian not the same.'
        raise Exception(msg)        
    # -------------------------------------------
    
    return sqrt(x)*exp(-x/T)*(2./sqrt(PI['value'])/sqrt(T)/T)


def conversion_to_SIUnits(values,units):
    """
    Converts array values from units to SI units.

    Input parameters:
    ---------------------------------------------
    values : array of values to be converted
    units  : array of units associated to values
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    SIvalues : array of converted values
    ---------------------------------------------
    """

    if isscalar(values) and isscalar(units):
        conversion_factor = get_conversion_factorSI(units)        
        SIvalues = values*conversion_factor         

    else:
        SIvalues_list = []

        if isscalar(units) and (isscalar(values) == False):
            conversion_factor = get_conversion_factorSI(units) 
            for value in values:  
                SIvalue = value*conversion_factor
                SIvalues_list.append(SIvalue)

        else:  
            for (value,unit) in zip(values,units):
                conversion_factor = get_conversion_factorSI(unit)        
                SIvalue = value*conversion_factor            
                SIvalues_list.append(SIvalue)

        SIvalues = array(SIvalues_list)

    return SIvalues


def get_conversion_factorSI(unit):
    """
    Gets conversion factor for a certain unit to SI units.

    Input parameters:
    ---------------------------------------------
    unit  : unit associated to be converted
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    conversion_factor : conversion factor to SI units based on unit
    ---------------------------------------------
    """

    conversion_factors = {'MeV': 1.60217653e-13, 'keV': 1.60217653e-16, 'eV': 1.60217653e-19, \
                          'ns': 1.0e-9, \
                          'kg': 1.0, \
                          'm': 1.0, 'dm': 0.1, 'cm': 0.01, 'mm': 0.001,\
                          'b': 1.0e-28, 'mb': 1.0e-31}

    if conversion_factors.get(unit) == None:
        msg =  "Error: Conversion factor", unit, "not implemented."
        raise Exception(msg)
    else:
        conversion_factor = conversion_factors[unit]

    return conversion_factor



def conversion_to_newUnits(values,unitsnew,unitsold):
    """
    Converts array values from units to SI units.

    Input parameters:
    ---------------------------------------------
    values   : array of values to be converted
    unitsold : array of units associated to values
    unitsnew : array of units, values should be
            converted to
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    convertedvalues : array of converted values
    ---------------------------------------------
    """

    if len(shape(values))<1 and isscalar(unitsold):
        conversion_factor = get_conversion_factor(unitsnew,unitsold)        
        convertedvalues = values*conversion_factor         

    else:
        convertedvalues_list = []

        if isscalar(unitsold) and (isscalar(values) == False):
            conversion_factor = get_conversion_factor(unitsnew,unitsold) 
            for value in values:  
                convertedvalue = value*conversion_factor
                convertedvalues_list.append(convertedvalue)

        else:  
            for (value,unitold) in zip(values,unitsold):
                conversion_factor = get_conversion_factor(unitsnew,unitold)        
                convertedvalue = value*conversion_factor            
                convertedvalues_list.append(convertedvalue)

        convertedvalues = array(convertedvalues_list)

    return convertedvalues



def get_conversion_factor(unitnew,unitold):
    """
    Gets conversion factor for a certain unit to another.

    Input parameters:
    ---------------------------------------------
    unitnew : unit to be converted to
    unitold : unit to be converted
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    conversion_factor : conversion factor
    ---------------------------------------------
    """

    conversion_factors = {'eV to MeV': 1.0e-6,\
                         'MeV to eV': 1.0e6}

    unitid = unitold + " to " + unitnew
    if conversion_factors.get(unitid) == None:
        msg =  "Error: Conversion factor", unitid, "not implemented."
        raise Exception(msg)
    else:
        conversion_factor = conversion_factors[unitid]

    return conversion_factor




def derivativefirst_linlininterpolation(lattice,data):
    """
    Calculates first derivative of lattice and
    data assuming lin-lin-interpolation of data.

    Input parameters:
    ---------------------------------------------
    lattice : array of lattice, must be larger 
         than 1.
    data   : array of data, where derivative is
        executed, size must be larger than 1.
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    lattice_derivative: lattice of derivative.
        Dimension is 1 entry smaller.
    data_derivative: first derivative of data
        Dimension is 1 entry smaller.
    ---------------------------------------------
    """

    # error messages if wrong lattice size ------
    if shape(lattice)[0] < 2 or shape(data)[0] < 2:
        raise Exception("Error in derivativefirst_linlininterpolation: lattice size must be at least 2.")

    if shape(lattice)[0] != shape(data)[0]:
        raise Exception("Error in derivativefirst_linlininterpolation: lattice and data must have the same shape.")
    # -------------------------------------------

    # initializing array ------------------------
    dim = shape(lattice)[0]

    lattice_derivative = zeros(dim-1,dtype='float')
    data_derivative    = zeros(dim-1,dtype='float')
    # -------------------------------------------

    # getting numerical derivative --------------
    if shape(lattice)[0] == 2:
        lattice_derivative =  (lattice[1]+lattice[0])*0.5
    else:
        lattice_derivative = (lattice[1:]+lattice[:-1])*0.5

    for index in arange(0,dim-1):
        data_derivative[index] = (data[index+1]-data[index])/(lattice[1+index]-lattice[index])
        if (lattice[1+index]-lattice[index])==0:
            print(lattice[1+index])
    # -------------------------------------------        

    return lattice_derivative,data_derivative
