# This module contains functions to manage the reference data

import os
currentdir = os.getcwd()
BASE_DIR = currentdir[0:currentdir.find("src")]
print(BASE_DIR)
REFERENCE_LIBRARY_FILE = BASE_DIR + "Data/ReferenceData_dict.txt"
print("REFERENCE_LIBRARY_FILE: " + str(REFERENCE_LIBRARY_FILE))


#import _pickle as cPickle
from pickle import dump, load
from numpy import loadtxt, shape, zeros


def get_reference_datacor(ref_lib):
    """
    Identifies file for reference data and reads
    them from the file.

    Input parameters:
    ---------------------------------------------
    ref_lib : library for reference, must contain
          information about the isotope, quantity,
          reaction and an identifier to be 
          identified. 
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    dataunccor_lib : library with tags 'Lattice',
          'Data', 'RelUnc' and 'Cor' containing
          associated information.
    ---------------------------------------------
    """

    print(ref_lib)
    ref_filename = identify_reference(ref_lib)
    dataunccor_lib = read_ref_file(ref_lib,ref_filename)

    return dataunccor_lib

def read_ref_file(ref_lib,ref_filename):
    """
    Identifies file for reference data.

    Input parameters:
    ---------------------------------------------
    ref_lib : library for reference, must contain
          information about the quantity.
    ref_filename : path to data and correlation
          matrix file implemented.
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    dataunccor_lib : data for the reference.
    tag 'lattice' 
    tag 'lattice2'
    tag 'lattice_unit'
    tag 'data'
    tag 'data_unit' : unit of data
    tag 'rel_unc' 
    tag 'cor': filename for the cor
    tag 'latticeunc' : lattice for rel_unc and
       cor if different from lattice
    tag 'latticeunc_unit' : unit of latticeunc 
    ---------------------------------------------
    """
    from numpy import interp, shape, arange
    import BasicPhysicsFunctions as BPF
    from importlib import reload
    reload(BPF)
    
    dataunccor_lib = {}

    if ref_lib['quantity'] == 'PFNS' :
        # reading data ===========================
        data = loadtxt(ref_filename['data_file_name'])
        dataunccor_lib['lattice'] = data[:,0]
        dataunccor_lib['lattice2'] = data[:,1]  
        dataunccor_lib['data'] = data[:,2]        
        dataunccor_lib['rel_unc'] = data[:,3]
        # ========================================

        # reading units ==========================
        f = open(ref_filename['data_file_name'])
        lines = f.readlines()
        f.close()

        dataunccor_lib['lattice_unit'] = lines[0].split()[2]
        dataunccor_lib['data_unit'] = lines[0].split()[6]

        del lines
        # ========================================


        # reading correlation matrix =============
        f = open(ref_filename['cor_file_name'])
        lines = f.readlines()
        lines = lines[1:]
        f.close()

        dim = shape(data)[0]
        cor = zeros([dim,dim],dtype = float)

        index1 = 0
        index2 = 0
        for line in lines:
            if line.split():
                cor[index1,index2] = float(line.split()[2])
                if index2 == dim-1:
                    index1 = index1 + 1
                    index2 = 0
                else:
                    index2 = index2 +1
                    
        del lines
        
        dataunccor_lib['cor'] = cor
        # ========================================

    elif ref_lib['quantity'] == 'cs' or ref_lib['quantity'] == 'nubar-tot' or ref_lib['quantity'] == 'nubar-prompt':
        # reading data ===========================
        data = loadtxt(ref_filename['data_file_name'])
        if ref_lib['reaction'] == '0,f' or ref_lib['reaction'] == 'n_th,f':
            dataunccor_lib['lattice'] = data[0]
            dataunccor_lib['data'] = data[1]       
        else:
            dataunccor_lib['lattice'] = data[:,0]
            dataunccor_lib['data'] = data[:,1]   
        
        dim = shape(data)[0]
        dataunccor_lib['rel_unc'] = zeros(dim,dtype = float)
        # ========================================

        # reading units ==========================
        f = open(ref_filename['data_file_name'])
        lines = f.readlines()
        f.close()

        dataunccor_lib['lattice_unit'] = lines[0].split()[2]
        dataunccor_lib['data_unit'] = lines[0].split()[4]

        del lines
        # ========================================


        # reading correlation matrix =============
        cor = zeros([dim,dim],dtype = float)

        if ref_filename['cor_file_name'] != 'none':
            # reading relative uncertainties ------
            data = loadtxt(ref_filename['cor_file_name'][1],skiprows=4)
            
            if ref_lib['reaction'] == '0,f' or ref_lib['reaction'] == 'n_th,f':
                Eincrelunc = data[0]
                RelUnc  = data[1]
            else:
                Eincrelunc = data[:,0]
                RelUnc  = data[:,1]
            # -------------------------------------
            
            # reading units of rel. unc. ----------
            f = open(ref_filename['cor_file_name'][1])
            lines = f.readlines()
            f.close()

            Eincrelunc_unit = lines[3].split()[1]
            # -------------------------------------

            if ref_lib['reaction'] == '0,f' or ref_lib['reaction'] == 'n_th,f':
                latticecor1 = [0]
                latticecor2 = [0]
                cor = [1]
            else:
                
                # reading correlation matrix ----------
                data = loadtxt(ref_filename['cor_file_name'][0],skiprows=4)
            
                latticecor1 = data[0,:]
                latticecor2 = data[1,:]
                cor = data[2:,:]/100.0

		# enforcing symmetry
                for index1 in arange(0,shape(cor)[0]):
                    for index2 in arange(index1+1,shape(cor)[0]):
                        cor[index1,index2] = cor[index2,index1]
                # ---------------------------------------

                # reading units of cor einc -------------
                f = open(ref_filename['cor_file_name'][0])
                lines = f.readlines()
                f.close()

                Einccor_unit = lines[4].split()[2]
                # ---------------------------------------

                # checking if energy grids are the same--
                if sum(latticecor1-latticecor2) != 0 :
                    msg = 'Error in '+ref_filename['cor_file_name'][0]+'. Covariance lattice of x1 and x2 not the same.'
                    raise(msg)
                
                if Eincrelunc_unit != Einccor_unit :
                    Eincrelunc = BPF.conversion_to_newUnits(Eincrelunc,Einccor_unit,Eincrelunc_unit)
                                
                if sum(latticecor1-Eincrelunc) != 0 :
                    RelUnc = array(interp(latticecor1,Eincrelunc,RelUnc))

                if dim != max(shape(latticecor1)) or sum(latticecor1-dataunccor_lib['lattice']):
                    dataunccor_lib['latticeunc'] = latticecor1
                    dataunccor_lib['latticeunc_unit'] = Einccor_unit
            
            dataunccor_lib['rel_unc'] = RelUnc
        # ---------------------------------------


        dataunccor_lib['cor'] = cor
        # ========================================  
    else:
        msg = "Error: Cannot read files of the quantity: " + str(ref_lib['quantity'])
        raise Exception(msg)                  

    return dataunccor_lib


def identify_reference(ref_lib):
    """
    Identifies file for reference data.

    Input parameters:
    ---------------------------------------------
    ref_lib : library for reference, must contain
          information about the isotope, quantity,
          reaction and an identifier to be 
          identified. 
    ---------------------------------------------

    Output parameters:
    ---------------------------------------------
    ref_filename : file names of the reference.
    tag 'data': filename for the data
    tag 'cor': filename for the cor
    ---------------------------------------------
    """
    ref_filename = {}

    all_references = load(open(REFERENCE_LIBRARY_FILE,'rb'))

    # finding filenames =========================
    for one_ref in all_references:
        comp_index = 0
        for key in ref_lib.keys():
            if one_ref[key] == ref_lib[key]:
                comp_index += 1
                if comp_index >= 4:
                    ref_filename['data_file_name'] = one_ref['data_file']
                    ref_filename['cor_file_name'] = one_ref['cor_file']
                    break
    # =============================================

    # error message if no reference found =========
    if ref_filename == {}:
        msg = "Error: Reference data not found."
        raise Exception(msg)        
    # =============================================

    return ref_filename



def create_reference_library():
    """
    Creates reference library binary file in REFERENCE_LIBRARY_FILE

    To execute: 
    add reference dictionary below and then run in ipython-notebook
    import Manage_ReferenceData as MR
    from importlib import reload
    reload(MR)
    MR.create_reference_library()
    
    """

    # reference files ========================================================
    referencetest = {'isotope': 'test_S-NNN','quantity': 'test_quantity', 'reaction': 'test_reaction',\
                     'identifier':'test_name','data_file':'test_data.dat','cor_file':'test_cor.dat'}
    
    reference_b10_VII1_nfcs = {'isotope': 'B-10','quantity': 'cs', 'reaction': 'n,a','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/B-10/n-005-B-010_nacs_ENDFBVII1.dat',\
                               'cor_file':[BASE_DIR+'Data/B-10/n-005-B-010_nacsCor_ENDFBVII1.dat',\
                                           BASE_DIR+'Data/B-10/n-005-B-010_nacsRelUnc_ENDFBVII1.dat']}
    
    reference_Li6_VII1_nfcs = {'isotope': 'Li-6','quantity': 'cs', 'reaction': 'n,t','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/Li-6/n-003-Li-006_ntcs_ENDFBVII1.dat',\
                               'cor_file':[BASE_DIR+'Data/Li-6/n-003-Li-006_ntcsCor_ENDFBVII1.dat',\
                                           BASE_DIR+'Data/Li-6/n-003-Li-006_ntcsRelUnc_ENDFBVII1.dat']}
                          
    reference_pu238_VIII0_nubartot = {'isotope': 'Pu-238','quantity': 'nubar-tot', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-238/n-094-Pu-238_nubartot_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu238_VIII0_nubarprompt = {'isotope': 'Pu-238','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-238/n-094-Pu-238_nubarprompt_ENDFBVIII0.dat',\
                          'cor_file':'none'}
    
    reference_pu239_VII1_nfcs = {'isotope': 'Pu-239','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_nfcs_ENDFBVII1.dat',\
                          'cor_file':'none'}
    
    reference_pu239_VIII0_nfcs = {'isotope': 'Pu-239','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_nfcs_ENDFBVIII0.dat',\
                          'cor_file':'none'}
    
    reference_pu239_VIII1_nfcs = {'isotope': 'Pu-239','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_nfcs_ENDFBVIII1.txt',\
                          'cor_file':[BASE_DIR+'Data/Pu-239/n-094-Pu-239_nfcsCor_ENDFBVIII1.txt',\
                                      BASE_DIR+'Data/Pu-239/n-094-Pu-239_nfcsRelUnc_ENDFBVIII1.txt']}
    
    reference_pu240_VIII0_nfcs = {'isotope': 'Pu-240','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-240/n-094-Pu-240_nfcs_ENDFBVIII0.dat',\
                          'cor_file':'none'}
    
    reference_pu240_PARADIGM_nfcs = {'isotope': 'Pu-240','quantity': 'cs', 'reaction': 'n,f','identifier':'PARADIGM',\
                          'data_file':BASE_DIR+'Data/Pu-240/n-094-Pu-240_nfcs_PARADIGM.dat',\
                          'cor_file':'none'}
                          
    reference_pu239_VIII0_nubartot = {'isotope': 'Pu-239','quantity': 'nubar-tot', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_nubartot_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu239_VIII0_nubarprompt = {'isotope': 'Pu-239','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_nubarprompt_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu240_VIII0_nubartot = {'isotope': 'Pu-240','quantity': 'nubar-tot', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-240/n-094-Pu-240_nubartot_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu240_VIII0_nubarprompt = {'isotope': 'Pu-240','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-240/n-094-Pu-240_nubarprompt_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu241_VIII0_nubartot = {'isotope': 'Pu-241','quantity': 'nubar-tot', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-241/n-094-Pu-241_nubartot_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu241_VIII0_nubarprompt = {'isotope': 'Pu-241','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-241/n-094-Pu-241_nubarprompt_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu242_VIII0_nubartot = {'isotope': 'Pu-242','quantity': 'nubar-tot', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-242/n-094-Pu-242_nubartot_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_pu242_VIII0_nubarprompt = {'isotope': 'Pu-242','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Pu-242/n-094-Pu-242_nubarprompt_ENDFBVIII0.dat',\
                          'cor_file':'none'}
                          
    reference_cf252_std2018_nubartot = {'isotope': 'Cf-252','quantity': 'nubar-tot', 'reaction': '0,f','identifier':'Standard2018',\
                          'data_file':BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubartot_standard2018.dat',\
                                        'cor_file':[BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubartotCor_standard2018.dat',\
                                        BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubartotRelUnc_standard2018.dat']}
    
    reference_cf252_VIII0_nubarprompt = {'isotope': 'Cf-252','quantity': 'nubar-prompt', 'reaction': '0,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubarprompt_VIII0.dat',\
                                        'cor_file':[BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubarpromptCor_VIII0.dat',\
                                        BASE_DIR+'Data/Cf-252/n-098-Cf-252_nubarpromptRelUnc_VIII0.dat']}

    reference_u235_VIII0thermal_nubarprompt = {'isotope': 'U-235','quantity': 'nubar-prompt', 'reaction': 'n_th,f','identifier':'ENDF/B-VIII.0_thermal',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nubarpromptthermal_VIII0.dat',\
                                        'cor_file':[BASE_DIR+'Data/U-235/n-092-U-235_nubarpromptthermalCor_VIII0.dat',\
                                        BASE_DIR+'Data/U-235/n-092-U-235_nubarpromptthermalRelUnc_VIII0.dat']}

    reference_u235_std2018_nubartot = {'isotope': 'U-235','quantity': 'nubar-tot', 'reaction': 'n_th,f','identifier':'Standard2018',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nubartotthermal_VIII0.dat',\
                                        'cor_file':[BASE_DIR+'Data/U-235/n-092-U-235_nubartotthermalCor_VIII0.dat',\
                                        BASE_DIR+'Data/U-235/n-092-U-235_nubartotthermalRelUnc_VIII0.dat']}
    
    reference_u232_VIII0_nubarprompt = {'isotope': 'U-232','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/U-232/n-092-U-232_nubarprompt_ENDFBVIII1.dat',\
                                        'cor_file':'none'}
    
    reference_u233_VIII0_nubarprompt = {'isotope': 'U-233','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/U-233/n-092-U-233_nubarprompt_ENDFBVIII0.dat',\
                                        'cor_file':'none'}
    
    reference_u234_VIII0_nubarprompt = {'isotope': 'U-234','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/U-234/n-092-U-234_nubarprompt_ENDFBVIII1.dat',\
                                        'cor_file':'none'}
    
    reference_u235_VIII0_nubarprompt = {'isotope': 'U-235','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nubarprompt_ENDFBVIII0.dat',\
                                        'cor_file':'none'}
    
    reference_u236_VIII0_nubarprompt = {'isotope': 'U-236','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/U-236/n-092-U-236_nubarprompt_ENDFBVIII1.dat',\
                                        'cor_file':'none'}
    
    reference_u238_VIII0_nubarprompt = {'isotope': 'U-238','quantity': 'nubar-prompt', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/U-238/n-092-U-238_nubarprompt_ENDFBVIII0.dat',\
                                        'cor_file':'none'}
    
    reference_u235_VII1_nfcs = {'isotope': 'U-235','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nfcs_ENDFBVII1.dat',\
                                'cor_file':[BASE_DIR+'Data/U-235/n-092-U-235_nfcsCor_ENDFBVII1.dat',\
                                            BASE_DIR+'Data/U-235/n-092-U-235_nfcsRelUnc_ENDFBVII1.dat']}

    reference_u235_VIII0_nfcs = {'isotope': 'U-235','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nfcs_ENDFBVIII0.dat',\
                                'cor_file':[BASE_DIR+'Data/U-235/n-092-U-235_nfcsCor_ENDFBVIII0.dat',\
                                            BASE_DIR+'Data/U-235/n-092-U-235_nfcsRelUnc_ENDFBVIII0.dat']}

    reference_u235_VIII1_nfcs = {'isotope': 'U-235','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/U-235/n-092-U-235_nfcs_ENDFBVIII1.txt',\
                                'cor_file':[BASE_DIR+'Data/U-235/n-092-U-235_nfcsCor_ENDFBVIII1.txt',\
                                            BASE_DIR+'Data/U-235/n-092-U-235_nfcsRelUnc_ENDFBVIII1.txt']}
    
    reference_u238_VII1_nfcs = {'isotope': 'U-238','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/U-238/n-092-U-238_nfcs_ENDFBVII1.dat',\
                                'cor_file':[BASE_DIR+'Data/U-238/n-092-U-238_nfcsCor_ENDFBVII1.dat',\
                                            BASE_DIR+'Data/U-238/n-092-U-238_nfcsRelUnc_ENDFBVII1.dat']}
    
    reference_u238_VIII0_nfcs = {'isotope': 'U-238','quantity': 'cs', 'reaction': 'n,f','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/U-238/n-092-U-238_nfcs_ENDFBVIII0.txt',\
                                'cor_file':'none'}
    
    reference_u238_VII1_ngcs = {'isotope': 'U-238','quantity': 'cs', 'reaction': 'n,g','identifier':'ENDF/B-VII.1',\
                          'data_file':BASE_DIR+'Data/U-238/n-092-U-238_ngcs_ENDFBVII1.dat',\
                                'cor_file':[BASE_DIR+'Data/U-238/n-092-U-238_ngcsCor_ENDFBVII1.dat',\
                                            BASE_DIR+'Data/U-238/n-092-U-238_ngcsRelUnc_ENDFBVII1.dat']}
    
    reference_u238_VIII1_ngcs = {'isotope': 'U-238','quantity': 'cs', 'reaction': 'n,g','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/U-238/n-092-U-238_ngcs_ENDFBVIII1.txt',\
                                'cor_file':'none'}
    
    reference_V51_VIII0_nfcs = {'isotope': 'V-51','quantity': 'cs', 'reaction': 'n,tot','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/V-51/n-023-V-051_ntotcs_ENDFBVIII.0.dat',\
                          'cor_file':'none'}
    
    reference_V51_CoH_nfcs = {'isotope': 'V-51','quantity': 'cs', 'reaction': 'n,tot','identifier':'CoH',\
                          'data_file':BASE_DIR+'Data/V-51/n-023-V-051_ntotcs_CoH.dat',\
                          'cor_file':'none'}
    
    reference_Cu63_VIII0_npcs = {'isotope': 'Cu-63','quantity': 'cs', 'reaction': 'n,p','identifier':'ENDF/B-VIII.0',\
                          'data_file':BASE_DIR+'Data/Cu-63/n-029-Cu-063_npcs_ENDFBVIII.0.dat',\
                          'cor_file':'none'}
    
    reference_Cu63_VIII1_ngcs = {'isotope': 'Cu-63','quantity': 'cs', 'reaction': 'n,gamma','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-63/n-029-Cu-063_ngcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cu63_VIII1_nelcs = {'isotope': 'Cu-63','quantity': 'cs', 'reaction': 'n,el','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-63/n-029-Cu-063_nelcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cu65_VIII1_ngcs = {'isotope': 'Cu-65','quantity': 'cs', 'reaction': 'n,gamma','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-65/n-029-Cu-065_ngcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cu65_VIII1_npcs = {'isotope': 'Cu-65','quantity': 'cs', 'reaction': 'n,p','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-65/n-029-Cu-065_npcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cu65_VIII1_nelcs = {'isotope': 'Cu-65','quantity': 'cs', 'reaction': 'n,el','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-65/n-029-Cu-065_nelcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cuel_VIII1_ngcs = {'isotope': 'Cu-el','quantity': 'cs', 'reaction': 'n,gamma','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-el/n-029-Cu-000_ngcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cuel_VIII1_ntotcs = {'isotope': 'Cu-el','quantity': 'cs', 'reaction': 'n,tot','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-el/n-029-Cu-000_ntotcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_Cuel_VIII1_nelcs = {'isotope': 'Cu-el','quantity': 'cs', 'reaction': 'n,el','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Cu-el/n-029-Cu-000_nelcs_ENDFBVIII1.txt',\
                          'cor_file':'none'}
    
    reference_pu239_VIII1_ngammacs = {'isotope': 'Pu-239','quantity': 'cs', 'reaction': 'n,gamma','identifier':'ENDF/B-VIII.1',\
                'data_file':BASE_DIR+'Data/Pu-239/n-094-Pu-239_ngammacs_ENDFBVIII1.dat',\
                'cor_file':'none'}
    
    reference_ti48_VIII1_ntotcs = {'isotope': 'Ti-48','quantity': 'cs', 'reaction': 'n,tot','identifier':'ENDF/B-VIII.1',\
                'data_file':BASE_DIR+'Data/Ti-48/n-022-Ti-048_ntotcs_ENDFBVIII1.txt',\
                'cor_file':'none'}
    
    reference_ti48_VIII1_ngcs = {'isotope': 'Ti-48','quantity': 'cs', 'reaction': 'n,g','identifier':'ENDF/B-VIII.1',\
                'data_file':BASE_DIR+'Data/Ti-48/n-022-Ti-048_ngcs_ENDFBVIII1.txt',\
                'cor_file':'none'}

    reference_au197_VII1_ngcs = {'isotope': 'Au-197','quantity': 'cs', 'reaction': 'n,gamma','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/Au-197/n-097-Au-197_ngcs_ENDFBVIII1.txt',\
                                'cor_file':[BASE_DIR+'Data/Au-197/n-097-Au-197_ngcsCor_ENDFBVIII1.txt',\
                                            BASE_DIR+'Data/Au-197/n-097-Au-197_ngcsRelUnc_ENDFBVIII1.txt']}

    reference_s32_VII1_npcs = {'isotope': 'S-32','quantity': 'cs', 'reaction': 'n,p','identifier':'ENDF/B-VIII.1',\
                          'data_file':BASE_DIR+'Data/S-32/n-016-S-032_npcs_ENDFBVIII1.txt',\
                                'cor_file':[BASE_DIR+'Data/S-32/n-016-S-032_npcsCor_ENDFBVIII1.txt',\
                                            BASE_DIR+'Data/S-32/n-016-S-032_npcsRelUnc_ENDFBVIII1.txt']}
    
    # =========================================================================

    # writing to file =========================================================
    reftotal =[referencetest,reference_pu239_VII1_nfcs,reference_u235_VII1_nfcs,reference_pu239_VIII0_nubartot,reference_pu239_VIII0_nubarprompt,\
               reference_cf252_std2018_nubartot,reference_cf252_VIII0_nubarprompt,reference_u235_VIII0thermal_nubarprompt,reference_u235_VIII0_nubarprompt,\
               reference_pu238_VIII0_nubartot,reference_pu238_VIII0_nubarprompt,reference_pu240_VIII0_nubartot,reference_pu240_VIII0_nubarprompt,\
               reference_pu241_VIII0_nubartot,reference_pu241_VIII0_nubarprompt,reference_pu242_VIII0_nubartot,reference_pu242_VIII0_nubarprompt,\
               reference_u235_std2018_nubartot,reference_pu239_VIII0_nfcs,reference_u235_VIII0_nfcs, reference_u238_VIII0_nubarprompt,reference_pu240_VIII0_nfcs,
               reference_V51_VIII0_nfcs,reference_b10_VII1_nfcs,reference_Li6_VII1_nfcs,reference_u238_VII1_nfcs,reference_u238_VIII0_nfcs,reference_u238_VII1_ngcs,reference_Cu63_VIII0_npcs,reference_V51_CoH_nfcs,reference_pu240_PARADIGM_nfcs,reference_u232_VIII0_nubarprompt,reference_u233_VIII0_nubarprompt,reference_u234_VIII0_nubarprompt,\
              reference_u236_VIII0_nubarprompt,reference_pu239_VIII1_nfcs,reference_pu239_VIII1_ngammacs,reference_ti48_VIII1_ntotcs,\
              reference_Cu63_VIII1_ngcs,reference_Cu65_VIII1_ngcs,reference_au197_VII1_ngcs,reference_u235_VIII1_nfcs,reference_ti48_VIII1_ngcs,reference_Cuel_VIII1_ngcs,reference_s32_VII1_npcs,reference_Cu65_VIII1_npcs,reference_Cuel_VIII1_ntotcs, reference_Cu63_VIII1_nelcs, reference_Cu65_VIII1_nelcs,\
               reference_Cuel_VIII1_nelcs,reference_u238_VIII1_ngcs]

    dump( reftotal, open(REFERENCE_LIBRARY_FILE, 'wb' ) )
    # =========================================================================
    return None

