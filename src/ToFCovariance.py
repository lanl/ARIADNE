from CorMatrixShapes import identify_corshape
from MatrixFunctions import checking_cov
from numpy import zeros, array, interp, isscalar, diag, ones, shape
import BasicPhysicsFunctions as BPF


# ======================================================================================================
#       N.W. TOF COVARIANCE FUNCTIONS
# ======================================================================================================

def nucdata_covariance_from_energy_covariance(deriv_nucdata, energy_covariance):
    J = diag(deriv_nucdata)
    nucdata_covariance = J @ energy_covariance @ J.T    
    return nucdata_covariance

def energy_covariance_from_length_unc(energy, energy_unit, 
                                      length, length_unit,
                                      length_unc, length_unc_unit):
    
    from MathsPhysicsConstants import NEUTRON_MASS
    from numpy import sqrt, array, diag
    from BasicPhysicsFunctions import get_conversion_factorSI

    if length_unc_unit == '%':
        length_unc = length_unc/100 * length
        length_unc_unit = length_unit

    length_unc_SI = BPF.conversion_to_SIUnits(length_unc, length_unc_unit)
    length_val_SI = BPF.conversion_to_SIUnits(length, length_unit) 
    energy_SI = BPF.conversion_to_SIUnits(energy,energy_unit)
    neutron_mass_SI = BPF.conversion_to_SIUnits(NEUTRON_MASS['value'],NEUTRON_MASS['unit'])   

    tof_SI = length_val_SI / sqrt( 2*energy_SI / neutron_mass_SI ) 
    dEi_dL = neutron_mass_SI/tof_SI**2 * length_val_SI
    Jac = array([dEi_dL]).T
    CovE = Jac @ diag([length_unc_SI**2]) @ Jac.T

    CovE = CovE / (get_conversion_factorSI(energy_unit)**2)

    return CovE


def energy_covariance_from_t0_unc(energy, energy_unit, 
                                length, length_unit,
                                t0_unc, t0_unc_unit):
    
    from MathsPhysicsConstants import NEUTRON_MASS
    from numpy import sqrt, array, diag
    from BasicPhysicsFunctions import get_conversion_factorSI

    if t0_unc_unit == '%':
        raise ValueError("% uncertainty in t0 not yet implemented, consider implications if t0 normalized to 0")

    length_val_SI = BPF.conversion_to_SIUnits(length, length_unit) 
    energy_SI = BPF.conversion_to_SIUnits(energy,energy_unit)
    t0_unc_SI = BPF.conversion_to_SIUnits(t0_unc, t0_unc_unit)
    neutron_mass_SI = BPF.conversion_to_SIUnits(NEUTRON_MASS['value'],NEUTRON_MASS['unit'])   

    tof_SI = length_val_SI / sqrt( 2*energy_SI / neutron_mass_SI ) 
    dEi_dt0 = energy_SI*-2/tof_SI
    Jac = array([dEi_dt0]).T
    CovE = Jac @ diag([t0_unc_SI**2]) @ Jac.T

    CovE = CovE / (get_conversion_factorSI(energy_unit)**2)

    return CovE


# ======================================================================================================
# energy-covariance functions that accept an ARIADNE CLASS, options for ratio and no ratio are the same for ALL CLASSES
# ======================================================================================================
from typing import Protocol, Any

class AriadneClass(Protocol):
    isotope : Any
    general_info : Any
    identifier_iso_deriv1 : Any
    einc_unit : Any
    values_unit : Any
    calculate_energy_resolution : Any
    einc : Any
    edges_for_chw_integrated_derivative : Any
    no_val : Any
    do_energy_unc : Any
    do_tof_length_unc : Any
    do_tof_t0_unc : Any
    values : Any

  
def calculate_cov_eunc_NORATIO(ariadne_class: AriadneClass):
    """Calculates cs covariance matrix due to energy uncertainty for data measured NOT in ratio.

    Args:
        ariadne_class (AriadneClass): any ARIADNE class instance

    Returns:
        ndarray: cs covariance matrix due to energy uncertainty
    """
    from ToFCovariance import energy_covariance_from_length_unc, energy_covariance_from_t0_unc, nucdata_covariance_from_energy_covariance
    from ExperimentalModelling import broaden_with_energy_resolution_vector
    from CorMatrixShapes import identify_corshape
    from Manage_ReferenceData import get_reference_datacor
    from MatrixFunctions import checking_cov
    from numpy import zeros, array, interp, diag, outer
    import BasicPhysicsFunctions as BPF
            
    # WIP program for multiple Tref and Tiso (different Einc)
    # WIP program for multiple Tref and Tiso (different Einc)
    if not any( (ariadne_class.do_energy_unc, ariadne_class.do_tof_length_unc, ariadne_class.do_tof_t0_unc) ):
        return zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    # reads in nuclear data and calculates numerical derivative -------------
    # reads data ............................................................
    nucdata_lib = {'isotope': ariadne_class.isotope,'quantity': ariadne_class.general_info['quantity'], \
                    'reaction':  ariadne_class.general_info['reaction'],'identifier': ariadne_class.identifier_iso_deriv1}
    dataunccor = get_reference_datacor(nucdata_lib)
    # .......................................................................

    # converts nuclear data to units of isoptope data .......................
    if dataunccor['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],ariadne_class.einc_unit,dataunccor['lattice_unit'])

    if dataunccor['data_unit'] != ariadne_class.values_unit:
        dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],ariadne_class.values_unit,dataunccor['data_unit'])
    # .......................................................................            

    # -----------------------------------------------------------------------
    # calculates derivative of experimental model WRT energy ................
    # check for energy resolution and bin width corrections .................
    model_data = dataunccor['data']
    model_lattice = dataunccor['lattice']

    if ariadne_class.calculate_energy_resolution:
        if len(shape(ariadne_class.einc)) >= 1:
            energy_resolution_fine_grid = interp(model_lattice, ariadne_class.einc, ariadne_class.energy_resolution)
        else:
            energy_resolution_fine_grid = ones(shape(model_lattice)[0])*ariadne_class.energy_resolution
        energy_resolution_broadened_model = broaden_with_energy_resolution_vector(
            model_lattice, model_data, 
            model_lattice, energy_resolution_fine_grid, 
        )
        model_data = energy_resolution_broadened_model

    if ariadne_class.edges_for_chw_integrated_derivative is None:
        [latticederiv,dataderiv] = BPF.derivativefirst_linlininterpolation(model_lattice, model_data)
        deriv_nucdata = array(interp(ariadne_class.einc,latticederiv,dataderiv))
    
    else:
        model_at_edges = interp(ariadne_class.edges_for_chw_integrated_derivative, model_lattice, model_data)
        # derivative defined by fundamental theorem of calculus - see derivation in tof_unc paper by Walton, N.A.W.         
        [_, deriv_nucdata] = BPF.derivativefirst_linlininterpolation(ariadne_class.edges_for_chw_integrated_derivative, model_at_edges)

    ariadne_class.model_data = model_data
    ariadne_class.model_lattice = model_lattice
    # .......................................................................            

    # get energy/time/length uncertainties into absolute energy covariance ..
    absolute_energy_covariance = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.do_energy_unc:
        if ariadne_class.enerr_unc_unit == '%':
            energy_unc = ariadne_class.enerr_unc/100 * ariadne_class.einc
            energy_unc_unit = ariadne_class.einc
        else:
            if ariadne_class.enerr_unc_unit != ariadne_class.einc_unit:
                ariadne_class.enerr_unc = BPF.conversion_to_newUnits(ariadne_class.enerr_unc, ariadne_class.einc_unit, ariadne_class.enerr_unc_unit)
                ariadne_class.enerr_unc_unit = ariadne_class.einc_unit
            energy_unc = ariadne_class.enerr_unc

        sub_cor = identify_corshape(ariadne_class.no_val, ariadne_class.enerr_unc_type, ariadne_class.enerr_unc_type_arg) 
        absolute_energy_covariance +=  sub_cor * outer(energy_unc,energy_unc)

    if ariadne_class.do_tof_length_unc:
        absolute_energy_covariance += energy_covariance_from_length_unc(ariadne_class.einc, ariadne_class.einc_unit, 
                                                                ariadne_class.tof_length_val, ariadne_class.tof_length_unit,
                                                                ariadne_class.tof_length_unc, ariadne_class.tof_length_unc_unit)

    if ariadne_class.do_tof_t0_unc:
        absolute_energy_covariance += energy_covariance_from_t0_unc(ariadne_class.einc, ariadne_class.einc_unit, 
                                                            ariadne_class.tof_length_val, ariadne_class.tof_length_unit,
                                                            ariadne_class.tof_t0_unc, ariadne_class.tof_t0_unc_unit)
    # .......................................................................            
    
    # propagate energy covariance to effective cov on experiment values .....
    nucdata_covariance = nucdata_covariance_from_energy_covariance(deriv_nucdata, absolute_energy_covariance)
    
    if len(shape(ariadne_class.einc)) >= 1:
        relative_nucdata_covariance = diag(1/ariadne_class.values) @ (nucdata_covariance*100*100) @ diag(1/ariadne_class.values)
    else:
        relative_nucdata_covariance = (1/ariadne_class.values) * (nucdata_covariance*100*100) * (1/ariadne_class.values)
    # .......................................................................            
    
    # -----------------------------------------------------------------------
    # FINAL testing .........................................................
    msg = "Testing cov_eunc covariance matrix:"
    print(msg)
    checking_cov(relative_nucdata_covariance)
    print("")
    # .......................................................................

    return relative_nucdata_covariance


def calculate_cov_eunc_RATIO(ariadne_class: AriadneClass):
    """Calculates cs covariance matrix due to energy uncertainty for data measured IN ratio.

    Args:
        ariadne_class (AriadneClass): any ARIADNE class instance

    Returns:
        ndarray: cs covariance matrix due to energy uncertainty
    """     
    from CorMatrixShapes import identify_corshape
    from ToFCovariance import nucdata_covariance_from_energy_covariance, energy_covariance_from_length_unc, energy_covariance_from_t0_unc
    from ExperimentalModelling import get_integrated_model, broaden_with_energy_resolution_vector
    from Manage_ReferenceData import get_reference_datacor
    from MatrixFunctions import checking_cov
    from numpy import zeros, array, interp, diag, any, diff, outer
    import BasicPhysicsFunctions as BPF
    from scipy.interpolate import interp1d
            
    # WIP program for multiple Tref and Tiso (different Einc)
    # WIP program for multiple Tref and Tiso (different Einc)

    # reads in nuclear data and calculates numerical derivative -------------
    # reads data ............................................................
    nucdata_lib = {'isotope': ariadne_class.isotope,'quantity': ariadne_class.general_info['quantity'], \
                    'reaction':  ariadne_class.general_info['reaction'],'identifier': ariadne_class.identifier_iso_deriv1}
    dataunccor = get_reference_datacor(nucdata_lib)
    dataunccor_ref = get_reference_datacor(ariadne_class.reference)
    # .......................................................................

    # converts nuclear data to units of isoptope data .......................
    if dataunccor['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],ariadne_class.einc_unit,dataunccor['lattice_unit'])
        dataunccor['lattice_unit'] = ariadne_class.einc_unit

    if dataunccor_ref['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor_ref['lattice'] = BPF.conversion_to_newUnits(dataunccor_ref['lattice'],ariadne_class.einc_unit,dataunccor_ref['lattice_unit'])
        dataunccor_ref['lattice_unit'] = ariadne_class.einc_unit

    if dataunccor['data_unit'] != dataunccor_ref['data_unit']:
        dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],dataunccor_ref['data_unit'],dataunccor['data_unit'])
        dataunccor['data_unit'] = dataunccor_ref['data_unit']
    # .......................................................................      

    # -----------------------------------------------------------------------
    # calculates derivative of experimental model WRT energy ................
    # check for energy resolution and bin width corrections .................
    model_data = dataunccor['data']
    model_lattice = dataunccor['lattice']
    ref_model_data = dataunccor_ref['data']
    ref_model_lattice = dataunccor_ref['lattice']

    if ariadne_class.calculate_energy_resolution:

        model_data = broaden_with_energy_resolution_vector(
                model_lattice, model_data, 
                model_lattice, interp(model_lattice, ariadne_class.einc, ariadne_class.energy_resolution), 
        )
        ref_model_data = broaden_with_energy_resolution_vector(
                ref_model_lattice, ref_model_data, 
                ref_model_lattice, interp(ref_model_lattice, ariadne_class.einc, ariadne_class.energy_resolution), 
        )


    if ariadne_class.edges_for_chw_integrated_derivative is None:
        ratio_data =  model_data / array(interp(model_lattice,ref_model_lattice,ref_model_data))
        [latticederiv,dataderiv] = BPF.derivativefirst_linlininterpolation(model_lattice,ratio_data)
        deriv_nucdata = array(interp(ariadne_class.einc,latticederiv,dataderiv))
        
    else:
        model_at_edges = interp(ariadne_class.edges_for_chw_integrated_derivative, model_lattice, model_data)
        ref_at_edges = interp(ariadne_class.edges_for_chw_integrated_derivative, ref_model_lattice, ref_model_data)

        [_, model_diff] = BPF.derivativefirst_linlininterpolation(ariadne_class.edges_for_chw_integrated_derivative, model_at_edges)
        [_, ref_diff] = BPF.derivativefirst_linlininterpolation(ariadne_class.edges_for_chw_integrated_derivative, ref_at_edges)

        int_model = get_integrated_model(ariadne_class.edges_for_chw_integrated_derivative, interp1d(model_lattice, model_data)) \
                    / diff(ariadne_class.edges_for_chw_integrated_derivative)
        int_ref = get_integrated_model(ariadne_class.edges_for_chw_integrated_derivative, interp1d(ref_model_lattice, ref_model_data)) \
                    / diff(ariadne_class.edges_for_chw_integrated_derivative)
        # quotient rule for derivatives - see derivation in tof_unc paper by Walton, N.A.W.         
        deriv_nucdata = (int_ref*model_diff - int_model*ref_diff) / int_ref**2

    ariadne_class.model_data = model_data
    ariadne_class.model_lattice = model_lattice
    ariadne_class.ref_model_data = ref_model_data
    ariadne_class.ref_model_lattice = ref_model_lattice
    # .......................................................................            

    # get energy/time/length uncertainties into absolute energy covariance ..
    absolute_energy_covariance = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.do_energy_unc:
        if ariadne_class.enerr_unc_unit == '%':
            energy_unc = ariadne_class.enerr_unc/100 * ariadne_class.einc
            energy_unc_unit = ariadne_class.einc
        else:
            if ariadne_class.enerr_unc_unit != ariadne_class.einc_unit:
                ariadne_class.enerr_unc = BPF.conversion_to_newUnits(ariadne_class.enerr_unc, ariadne_class.einc_unit, ariadne_class.enerr_unc_unit)
                ariadne_class.enerr_unc_unit = ariadne_class.einc_unit
            energy_unc = ariadne_class.enerr_unc

        sub_cor = identify_corshape(ariadne_class.no_val, ariadne_class.enerr_unc_type, ariadne_class.enerr_unc_type_arg) 
        absolute_energy_covariance +=  sub_cor * outer(energy_unc,energy_unc)
        
    if ariadne_class.do_tof_length_unc:
        absolute_energy_covariance += energy_covariance_from_length_unc(ariadne_class.einc, ariadne_class.einc_unit, 
                                                                ariadne_class.tof_length_val, ariadne_class.tof_length_unit,
                                                                ariadne_class.tof_length_unc, ariadne_class.tof_length_unc_unit)

    if ariadne_class.do_tof_t0_unc:
        absolute_energy_covariance += energy_covariance_from_t0_unc(ariadne_class.einc, ariadne_class.einc_unit, 
                                                            ariadne_class.tof_length_val, ariadne_class.tof_length_unit,
                                                            ariadne_class.tof_t0_unc, ariadne_class.tof_t0_unc_unit)
    # .......................................................................            
    
    # propagate energy covariance to effective cov on experiment values .....
    if any( (ariadne_class.do_energy_unc, ariadne_class.do_tof_length_unc, ariadne_class.do_tof_t0_unc) ):
        nucdata_covariance = nucdata_covariance_from_energy_covariance(deriv_nucdata, absolute_energy_covariance)
        relative_nucdata_covariance = diag(1/ariadne_class.values) @ (nucdata_covariance*100*100) @ diag(1/ariadne_class.values)
    else:
        relative_nucdata_covariance = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    # .......................................................................            

    # -----------------------------------------------------------------------

    # Test and return -------------------------------------------------------
    # FINAL testing .........................................................
    msg = "Testing cov_eunc covariance matrix:"
    print(msg)
    checking_cov(relative_nucdata_covariance)
    print("")
    # .......................................................................

    return relative_nucdata_covariance





# ======================================================================================================
#       LEGACY energy covariance FUNCTIONS FOR REPRODUCIBILITY/REFERENCE
# ======================================================================================================

def calculate_cov_eunc_NORATIO_legacy(ariadne_class):
    """
    Calculates Eunc covariance matrix relative to cs.
    
    output variable:
    cov_eunc : Eunc covariances matrix relative to cs.
    """        
    from CorMatrixShapes import identify_corshape
    from Manage_ReferenceData import get_reference_datacor
    from MatrixFunctions import checking_cov
    from numpy import zeros, isscalar, array, interp, outer, isinf, arange, shape
    import BasicPhysicsFunctions as BPF
    from importlib import reload
    reload(BPF)
    from MathsPhysicsConstants import NEUTRON_MASS

    cov_eunc = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
            
    # WIP program for multiple Tref and Tiso (different Einc)
    # WIP program for multiple Tref and Tiso (different Einc)
    
    # check if Eout unc. are given in per-cent and converts that to unit of
    # Eout for calculation --------------------------------------------------
    if ariadne_class.enerr_unc_unit == '%':
        ariadne_class.enerr_unc = ariadne_class.enerr_unc/100
    else:
        msg = 'Error: conversion of units not implemented for Eunc.'
        raise Exception(msg)  
        
    if ariadne_class.enrsl_unc_unit == '%':
        ariadne_class.enrsl_unc = ariadne_class.enrsl_unc/100
    else:
        msg = 'Error: conversion of units not implemented for Eunc.'
        raise Exception(msg)  
    # -----------------------------------------------------------------------


    # reads in nuclear data and calculates numerical derivative -------------
    # reads data ............................................................
    nucdata_lib = {'isotope': ariadne_class.isotope,'quantity': ariadne_class.general_info['quantity'], \
                    'reaction':  ariadne_class.general_info['reaction'],'identifier': ariadne_class.identifier_iso_deriv1}
    dataunccor = get_reference_datacor(nucdata_lib)
    # .......................................................................

    # converts nuclear data to units of isoptope data .......................
    if dataunccor['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],ariadne_class.einc_unit,dataunccor['lattice_unit'])
        dataunccor['lattice_unit'] = ariadne_class.einc_unit

    #if dataunccor['data_unit'] != ariadne_class.values_unit:
    #    dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],ariadne_class.values_unit,dataunccor['data_unit'])
    # .......................................................................            

    # calculates numerical first derivative 
    [latticederiv,dataderiv] = BPF.derivativefirst_linlininterpolation(dataunccor['lattice'],dataunccor['data'])
    
    # interpolates first derivative 
    deriv_nucdata = array(interp(ariadne_class.einc,latticederiv,dataderiv))
    # -----------------------------------------------------------------------
    
    for index in arange(0,shape(deriv_nucdata)[0]):
        if isinf(deriv_nucdata[index]):
            deriv_nucdata[index] = 0
    
    # calculates cov_Eunc due to enerr --------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.no_val == 1 and ariadne_class.enerr_unc != 0:
        sub_cov[0,0] = ariadne_class.enerr_unc*ariadne_class.enerr_unc*deriv_nucdata*deriv_nucdata*ariadne_class.einc*ariadne_class.einc/\
                        (ariadne_class.values*ariadne_class.values)
    else:          
        if any(ariadne_class.enerr_unc): # non-zero uncertainties
            sub_cor = identify_corshape(ariadne_class.no_val,ariadne_class.enerr_unc_type,ariadne_class.enerr_unc_type_arg)
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = ariadne_class.enerr_unc[index1]*ariadne_class.enerr_unc[index2]*\
                                            deriv_nucdata[index1]*deriv_nucdata[index2]*\
                                            sub_cor[index1,index2]*ariadne_class.einc[index1]*ariadne_class.einc[index2]/\
                                            (ariadne_class.values[index1]*ariadne_class.values[index2])

        # testing ...........................................................
        msg = "Testing enerr covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    # calculates cov_Eunc due to enrsl --------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.no_val == 1 and ariadne_class.enrsl_unc != 0:
        sub_cov[0,0] = ariadne_class.enrsl_unc*ariadne_class.enrsl_unc*deriv_nucdata*deriv_nucdata*ariadne_class.einc*ariadne_class.einc/\
                        (ariadne_class.values*ariadne_class.values)
    else:          
        if any(ariadne_class.enrsl_unc): # non-zero uncertainties
            sub_cor = identify_corshape(ariadne_class.no_val,ariadne_class.enrsl_unc_type,ariadne_class.enrsl_unc_type_arg)
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = ariadne_class.enrsl_unc[index1]*ariadne_class.enrsl_unc[index2]*\
                                            deriv_nucdata[index1]*deriv_nucdata[index2]*\
                                            sub_cor[index1,index2]*ariadne_class.einc[index1]*ariadne_class.einc[index2]/\
                                            (ariadne_class.values[index1]*ariadne_class.values[index2])

        # testing ...........................................................
        msg = "Testing enrsl covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    # calculates cov_Eunc due to trsl ---------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.trsl_val != 0: # non-zero uncertainties

        # convert to SI units ---------------------------------------------------
        trsl_val_SI = BPF.conversion_to_SIUnits(ariadne_class.trsl_val,ariadne_class.trsl_unit)
        tof_length_val_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_val,ariadne_class.tof_length_unit) 
        einc_SI = BPF.conversion_to_SIUnits(ariadne_class.einc,ariadne_class.einc_unit)
        neutron_mass_SI = BPF.conversion_to_SIUnits(NEUTRON_MASS['value'],NEUTRON_MASS['unit'])   

        # convert derivative to SI units 

        lattice_SI = BPF.conversion_to_SIUnits(dataunccor['lattice'],dataunccor['lattice_unit'])
        data_SI = BPF.conversion_to_SIUnits(dataunccor['data'],dataunccor['data_unit'])
        [latticederiv_SI,dataderiv_SI] = BPF.derivativefirst_linlininterpolation(lattice_SI,data_SI)
        deriv_nucdata_SI = array(interp(einc_SI,latticederiv_SI,dataderiv_SI))
        shapecs_convserionfactor = BPF.conversion_to_SIUnits(1.0,dataunccor['data_unit'])
        # -----------------------------------------------------------------------
        
        # calculate cov ---------------------------------------------------------
        if isscalar(trsl_val_SI) and isscalar(tof_length_val_SI):
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = 8*pow(einc_SI[index1],1.5)*pow(einc_SI[index2],1.5)*\
                                                deriv_nucdata_SI[index1]*deriv_nucdata_SI[index2]*\
                                                pow(trsl_val_SI,2.0)/(neutron_mass_SI*pow(tof_length_val_SI,2.0)*\
                                                                    ariadne_class.values[index1]*ariadne_class.values[index2]*pow(shapecs_convserionfactor,2.0))
        else:
            msg = 'Error: non-scalar trsl, tof-length not implemented.'
            raise Exception(msg)            
            # -----------------------------------------------------------------------

        # testing ...........................................................
        msg = "Testing trsl covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    
    # calculates cov_Eunc due to tof_length ---------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.tof_length_unc != 0: # non-zero uncertainties
        
        # convert to SI units ---------------------------------------------------
        tof_length_unc_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_unc,ariadne_class.tof_length_unc_unit)
        tof_length_val_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_val,ariadne_class.tof_length_unit) 
        einc_SI = BPF.conversion_to_SIUnits(ariadne_class.einc,ariadne_class.einc_unit)

        # convert derivative to SI units 
        lattice_SI = BPF.conversion_to_SIUnits(dataunccor['lattice'],dataunccor['lattice_unit'])
        data_SI = BPF.conversion_to_SIUnits(dataunccor['data'],dataunccor['data_unit'])
        [latticederiv_SI,dataderiv_SI] = BPF.derivativefirst_linlininterpolation(lattice_SI,data_SI)
        deriv_nucdata_SI = array(interp(einc_SI,latticederiv_SI,dataderiv_SI))
        shapecs_convserionfactor = BPF.conversion_to_SIUnits(1.0,dataunccor['data_unit'])
        # -----------------------------------------------------------------------
        
        # calculate cov ---------------------------------------------------------
        if isscalar(tof_length_unc_SI) and isscalar(tof_length_val_SI):
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = 4*deriv_nucdata_SI[index1]*deriv_nucdata_SI[index2]*einc_SI[index1]*\
                                                einc_SI[index2]*pow(tof_length_unc_SI,2.0)/(pow(tof_length_val_SI,2.0)*\
                                                ariadne_class.values[index1]*ariadne_class.values[index2]*pow(shapecs_convserionfactor,2.0))
        else:
            msg = 'Error: non-scalar tof-length or tof-length unc. not implemented.'
            raise Exception(msg)            
        # -----------------------------------------------------------------------

        # testing ...........................................................
        msg = "Testing tof_length_val covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------

    
    # testing ---------------------------------------------------------------
    msg = "Testing cov_eunc covariance matrix:"
    print(msg)
    checking_cov(cov_eunc)
    print("")
    # -----------------------------------------------------------------------

    # freeing variables -----------------------------------------------------
    del sub_cov, sub_cor, deriv_nucdata,latticederiv,dataderiv, dataunccor
    del nucdata_lib 
    # -----------------------------------------------------------------------

    return cov_eunc*(100*100)


def calculate_cov_eunc_RATIO_legacy(ariadne_class):
    """
    Calculates Eunc covariance matrix relative to cs.
    
    output variable:
    cov_eunc : Eunc covariances matrix relative to cs.
    """        
    from CorMatrixShapes import identify_corshape
    from Manage_ReferenceData import get_reference_datacor
    from MatrixFunctions import checking_cov
    from numpy import zeros, isscalar, array, interp, diag,sqrt
    import BasicPhysicsFunctions as BPF
    from importlib import reload
    reload(BPF)
    from MathsPhysicsConstants import NEUTRON_MASS

    cov_eunc = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
            
    # WIP program for multiple Tref and Tiso (different Einc)
    # WIP program for multiple Tref and Tiso (different Einc)
    
    # check if Eout unc. are given in per-cent and converts that to unit of
    # Eout for calculation --------------------------------------------------
    if ariadne_class.enerr_unc_unit == '%':
        ariadne_class.enerr_unc = ariadne_class.enerr_unc/100
    else:
        msg = 'Error: conversion of units not implemented for Eunc.'
        raise Exception(msg)  
        
    if ariadne_class.enrsl_unc_unit == '%':
        ariadne_class.enrsl_unc = ariadne_class.enrsl_unc/100
    else:
        msg = 'Error: conversion of units not implemented for Eunc.'
        raise Exception(msg)  
    # -----------------------------------------------------------------------


    # reads in nuclear data and calculates numerical derivative -------------
    # reads data ............................................................
    nucdata_lib = {'isotope': ariadne_class.isotope,'quantity': ariadne_class.general_info['quantity'], \
                    'reaction':  ariadne_class.general_info['reaction'],'identifier': ariadne_class.identifier_iso_deriv1}
    dataunccor = get_reference_datacor(nucdata_lib)
    dataunccor_ref = get_reference_datacor(ariadne_class.reference)
    # .......................................................................

    # converts nuclear data to units of isoptope data .......................
    if dataunccor['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],ariadne_class.einc_unit,dataunccor['lattice_unit'])
        dataunccor['lattice_unit'] = ariadne_class.einc_unit

    if dataunccor_ref['lattice_unit'] != ariadne_class.einc_unit:
        dataunccor_ref['lattice'] = BPF.conversion_to_newUnits(dataunccor_ref['lattice'],ariadne_class.einc_unit,dataunccor_ref['lattice_unit'])
        dataunccor_ref['lattice_unit'] = ariadne_class.einc_unit

    if dataunccor['data_unit'] != dataunccor_ref['data_unit']:
        dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],dataunccor_ref['data_unit'],dataunccor['data_unit'])
        dataunccor['data_unit'] = dataunccor_ref['data_unit']
    # .......................................................................            

    # calculates ratio data 
    data_ref_intp = array(interp(dataunccor['lattice'],dataunccor_ref['lattice'],dataunccor_ref['data']))
    ratio_data =  dataunccor['data']/data_ref_intp

    # calculates numerical first derivative 
    [latticederiv,dataderiv] = BPF.derivativefirst_linlininterpolation(dataunccor['lattice'],ratio_data)
    
    # interpolates first derivative 
    deriv_nucdata = array(interp(ariadne_class.einc,latticederiv,dataderiv))
    # -----------------------------------------------------------------------
    
    # calculates cov_Eunc due to enerr --------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.no_val == 1 and ariadne_class.enerr_unc != 0:
        sub_cov[0,0] = ariadne_class.enerr_unc*ariadne_class.enerr_unc*deriv_nucdata*deriv_nucdata*ariadne_class.einc*ariadne_class.einc/\
                        (ariadne_class.values*ariadne_class.values)
    else:          
        if any(ariadne_class.enerr_unc): # non-zero uncertainties
            sub_cor = identify_corshape(ariadne_class.no_val,ariadne_class.enerr_unc_type,ariadne_class.enerr_unc_type_arg)
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = ariadne_class.enerr_unc[index1]*ariadne_class.enerr_unc[index2]*\
                                            deriv_nucdata[index1]*deriv_nucdata[index2]*\
                                            sub_cor[index1,index2]*ariadne_class.einc[index1]*ariadne_class.einc[index2]/\
                                            (ariadne_class.values[index1]*ariadne_class.values[index2])

        # testing ...........................................................
        msg = "Testing enerr covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    # calculates cov_Eunc due to enrsl --------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.no_val == 1 and ariadne_class.enrsl_unc != 0:
        sub_cov[0,0] = ariadne_class.enrsl_unc*ariadne_class.enrsl_unc*deriv_nucdata*deriv_nucdata*ariadne_class.einc*ariadne_class.einc/\
                        (ariadne_class.values*ariadne_class.values)
    else:          
        if any(ariadne_class.enrsl_unc): # non-zero uncertainties
            sub_cor = identify_corshape(ariadne_class.no_val,ariadne_class.enrsl_unc_type,ariadne_class.enrsl_unc_type_arg)
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = ariadne_class.enrsl_unc[index1]*ariadne_class.enrsl_unc[index2]*\
                                            deriv_nucdata[index1]*deriv_nucdata[index2]*\
                                            sub_cor[index1,index2]*ariadne_class.einc[index1]*ariadne_class.einc[index2]/\
                                            (ariadne_class.values[index1]*ariadne_class.values[index2])
        # testing ...........................................................
        msg = "Testing enrsl covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    # calculates cov_Eunc due to trsl ---------------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.trsl_val != 0: # non-zero uncertainties

        # convert to SI units ---------------------------------------------------
        trsl_val_SI = BPF.conversion_to_SIUnits(ariadne_class.trsl_val,ariadne_class.trsl_unit)
        tof_length_val_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_val,ariadne_class.tof_length_unit) 
        neutron_mass_SI = BPF.conversion_to_SIUnits(NEUTRON_MASS['value'],NEUTRON_MASS['unit'])
        conversionfactor_E = BPF.conversion_to_SIUnits(1.0,ariadne_class.einc_unit)
        
        # calculate cov ---------------------------------------------------------
        if isscalar(trsl_val_SI) and isscalar(tof_length_val_SI):
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = 8.0*pow(ariadne_class.einc[index1],1.5)*pow(ariadne_class.einc[index2],1.5)*\
                                                deriv_nucdata[index1]*deriv_nucdata[index2]*conversionfactor_E*\
                                                pow(trsl_val_SI,2.0)/(neutron_mass_SI*pow(tof_length_val_SI,2.0)*\
                                                                    ariadne_class.values[index1]*ariadne_class.values[index2])
        else:
            msg = 'Error: non-scalar trsl, tof-length not implemented.'
            raise Exception(msg)            
            # -----------------------------------------------------------------------

        # testing ...........................................................
        msg = "Testing trsl covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    
    # calculates cov_Eunc due to tof_length ---------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    if ariadne_class.tof_length_unc != 0: # non-zero uncertainties
        # convert to SI units ---------------------------------------------------
        tof_length_unc_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_unc,ariadne_class.tof_length_unc_unit)
        tof_length_val_SI = BPF.conversion_to_SIUnits(ariadne_class.tof_length_val,ariadne_class.tof_length_unit) 
        # -----------------------------------------------------------------------
        
        # calculate cov ---------------------------------------------------------
        if isscalar(tof_length_unc_SI) and isscalar(tof_length_val_SI):
            for index1 in range(0,ariadne_class.no_val):
                for index2 in range(0,ariadne_class.no_val):
                    sub_cov[index1,index2] = 4.0*deriv_nucdata[index1]*deriv_nucdata[index2]*ariadne_class.einc[index1]*\
                                                ariadne_class.einc[index2]*pow(tof_length_unc_SI,2.0)/(pow(tof_length_val_SI,2.0)*\
                                                ariadne_class.values[index1]*ariadne_class.values[index2])
        else:
            msg = 'Error: non-scalar trsl, tof-length not implemented.'
            raise Exception(msg)            
            # -----------------------------------------------------------------------

        # testing ...........................................................
        msg = "Testing trsl covariance matrix:"
        print(msg)
        checking_cov(sub_cov)
        print("")
        # ...................................................................
        
    cov_eunc += sub_cov
    # -----------------------------------------------------------------------
    
    
    # calculates cov_Eunc due to tof_length ---------------------------------
    sub_cov  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)
    sub_cor  = zeros([ariadne_class.no_val,ariadne_class.no_val],dtype = float)

    
    # testing ---------------------------------------------------------------
    msg = "Testing cov_eunc covariance matrix:"
    print(msg)
    checking_cov(cov_eunc)
    print("")
    # -----------------------------------------------------------------------

    # freeing variables -----------------------------------------------------
    del sub_cov, sub_cor, deriv_nucdata,latticederiv,dataderiv, dataunccor
    del nucdata_lib 
    # -----------------------------------------------------------------------

    return cov_eunc*(100*100)
