# This module contains functions to write out data
import os
currentdir = os.getcwd()
BASE_DIR = currentdir[0:currentdir.find("src")]

def write_partialunc_output(path_name,data):
    """
    Partial uncertainties are saved in a simple text-file.

    ------------------------------------------------------
    Input:
    path_name : string containing folder, where data should
    be stored.
    data : dictionary containing data array and header for
    the uncertainties
    ------------------------------------------------------

    ------------------------------------------------------
    Output:
    msg : stating path to file
    ------------------------------------------------------    
    """
    from numpy import shape

    
    # writing to file ------------------------------------
    filename = path_name + 'Partial_Unc.dat'
    f = open(filename,'w')

    print('#'+ ' '.join(map(str,data['header'])),file=f)
    
    if shape(shape(data['values']))[0] != 1:
        for index1 in range(0,shape(data['values'])[1]):
            print(' '.join(map(str,data['values'][:,index1])),file=f)
    else:
        print(' '.join(map(str,data['values'][:])),file=f)

    f.close()
    # ----------------------------------------------------

    msg = "Output file saved in file: " + str(filename)   
    print("")
     
    return msg

def write_cortofile(filepath_name,lattice,cor):
    """
    Correlation matrix is written to file in a way that it
    can be stored by gnuplot.

    ------------------------------------------------------
    Input:
    filepath_name : string containing file and folder name, 
    where data should be stored.
    lattice : lattices for cor
    cor : 2dim array containing correlation matrix
    ------------------------------------------------------

    ------------------------------------------------------
    Output:
    msg : stating path to file
    ------------------------------------------------------    
    """
    from numpy import shape

    # checking of input data -----------------------------
    if shape(shape(lattice))[0] != 0:
        if len(shape(cor)) != 2:
            msg = 'Error: Expect two dimensional array. But shape is: ' + str(shape(cor))
            raise Exception(msg)
        if shape(cor)[0] != shape(cor)[1]:
            msg = 'Error: Matrix must be square.'
            raise Exception(msg)
        if shape(cor)[0] != shape(lattice)[0]:
            msg = 'Error: Lattice and correlation matrix must have the same number of rows.'
            raise Exception(msg)        
    # ----------------------------------------------------

    
    # writing to file ------------------------------------
        f = open(filepath_name,'w')

        for index1 in range(0,shape(cor)[1]):
            for index2 in range(0,shape(cor)[1]):        
                print(lattice[index1],lattice[index2],cor[index1,index2],file=f)
                print(" ",file=f)
                
        f.close()
    # ----------------------------------------------------

    
    # writing to file easy to read for loaadtxt ------------------------------------
        f = open(filepath_name+"loadtxt",'w')
        
        print("# line 2: x-coordinate lattice, line 3: y-coordinate lattice starting line 4: correlation entries, each line a row",file=f)
        print(' '.join(map(str,lattice)),file=f)
        print(' '.join(map(str,lattice)),file=f)
        for index1 in range(0,shape(cor)[1]):
            print(" ".join(map(str,cor[index1,:])),file=f)
                
        f.close()
    # ----------------------------------------------------

        msg = "Output file saved in file: " + str(filepath_name) 
        print("")  
    else:
         msg = "No correlation matrix saved as only one data point given."       
     
    return msg

def write_xml_output(general_info,data):
    """
    Writes xml_output.
    
    Input:
    general_info : dictionary containing name, isotope, quantity, reaction
    and name of the output file.
    data : dictionary containing einc_unit, eout_unit, values_unit, no_einc,
    einc, eout, values, rel_unc and cor.
    
    output:
    msg : message stating path to file.
    """
    from lxml import etree
    from numpy import array, shape, ravel
    import datetime
    
    filename = general_info['output_file']
    
    # main directory ==================================================
    main = etree.Element('reactionSuite', \
                         projectile = general_info['reaction'][0], \
                         target = general_info['isotope'], \
                         formatVersion = 'GND 1.7', \
                         library = 'Estimated experimental data and covariances', \
                         version = "0.1", \
                         projectileFrame="lab" )
    # ==================================================================
    
    # styles entry =====================================================
    styles = etree.SubElement(main, 'styles')
    expunc = etree.SubElement(styles, 'exp_unc', \
                              name = "estimated experimental uncertainty",\
                              date = str(datetime.datetime.now().date()) )
    # ==================================================================
    
    # WIP: each experiment is a seperate xml-file, cross-covariances will
    # also be a separate file
    
    # documentation entry ==============================================
    documentations = etree.SubElement(main, 'documentations')
    doc_name = str(general_info['name'])+'_'+str(general_info['isotope'])\
    +'_'+str(general_info['quantity'])
    
    doc_exp_unc = etree.SubElement(documentations, 'documentation', \
                                   name = doc_name)   
    
    doc_exp_unc.text = str(general_info['documentation'])
    # ==================================================================
    
    # write data and covariances =======================================
    if general_info['quantity'] == 'PFNS':
        reaction = etree.SubElement(main, 'reaction', label=str(45), \
                outputChannel="n[emissionMode:'prompt']", ENDF_MT=str(18), \
                fissionGenre="total")
        outputChannel = etree.SubElement(reaction,'outputChannel',genre="NBody")
        product = etree.SubElement(outputChannel, 'product', name="n", label="n",\
                                   emissionMode="prompt")
        distributions = etree.SubElement(product, 'distributions')
        uncorrelated = etree.SubElement(distributions, 'uncorrelated',\
                                        style="exp_unc",productFrame="lab")
        multiD_XYs = etree.SubElement(uncorrelated,"multiD_XYs",\
                                      dimension="2",label="energy")
        
        # axis ----------------------------------------------------------
        axes = etree.SubElement(multiD_XYs,"axes")
        einc_axis = etree.SubElement(axes,'axis',index="2",label="energy_in",\
                                unit=data['einc_unit'] )
        eout_axis = etree.SubElement(axes,'axis',index="1",label="energy_out",\
                                unit=data['eout_unit'] )
        PFNS_axis = etree.SubElement(axes,'axis',index="0",label="P(energy_out|energy_in)",\
                                unit=data['values_unit'] )
        # ----------------------------------------------------------------
        
        # xys ------------------------------------------------------------
        indexEinc = 0
        for index2 in range(0,data['no_einc']):
            xys = etree.SubElement(multiD_XYs, "XYs",index=str(index2),\
                                   value=str(data['einc'][indexEinc]))
            if data['no_einc'] == 1:
                noval = shape(data['data'])[0]
            else:
                noval = 0
                while data['einc'][indexEinc] == data['einc'][indexEinc+noval]:
                    noval += 1
                    if (noval+indexEinc) == shape(data['data'])[0]:
                        break
            
            # values .....................................................
            val = etree.SubElement(xys, 'values', length = str(2*noval))
            val.text = ' '.join(map(str,ravel(list(zip(data['eout'][indexEinc:indexEinc+noval], \
                        data['data'][indexEinc:indexEinc+noval])))))
            # ............................................................
            
            # relative uncertainties .....................................
            uncs = etree.SubElement(xys, 'uncertainties')
            unc = etree.SubElement(uncs,'uncertainty', type = "variance",\
                                   relation = 'relative')
            val = etree.SubElement(unc, 'values', length = str(2*noval))
            val.text = ' '.join(map(str,ravel(list(zip(data['eout'][indexEinc:indexEinc+noval], \
                        data['rel_unc'][indexEinc:indexEinc+noval])))) )           
            # ............................................................                     
            
            indexEinc = indexEinc + noval           
        
        # ----------------------------------------------------------------
        
        # covariances ----------------------------------------------------
        cov = etree.SubElement(multiD_XYs, "covarianceMatrix", type = "cor")
        covmat = etree.SubElement(cov, "matrix", rows = str(shape(data['data'])[0]), \
                    columns = str(shape(data['data'])[0]), form = "symmetric")
        
        cor_data = data['cor']
        covflat = []
        for index1 in range(0,shape(data['data'])[0]):
            covflat = covflat + list(map(str,cor_data[index1,index1:]))
        
        covmat.text = ' '.join(covflat)
        #covmat.text = ' '.join(['1.0 ' '1.0'])  
        # ----------------------------------------------------------------
    elif general_info['quantity'] == 'cs' or general_info['quantity'] == 'nubar-tot' or general_info['quantity'] == 'nubar-prompt' : 
        msg = "WIP: not yet included. No Xml file produced."
        print(msg)
    else:
        msg = "Error: write xml  output not implemented for quantity" + str(general_info['quantity'])
        raise Exception(msg)                                      
    
    # ==================================================================
    
    
    # write xml file ===================================================
    xml_file = open(filename,'wb')
    main_file = etree.ElementTree(main)
    main_file.write(xml_file)
    # ==================================================================
    
    msg = "Output file saved in file: " + str(filename) 
    print("")       

    return msg

def write_json_output(general_info,data,features):
    """
    Writes json_output.
    
    Input:
    general_info : dictionary containing name, isotope, quantity, reaction
    and name of the output file.
    data : dictionary containing einc_unit, eout_unit, values_unit, no_einc,
    einc, eout, values, rel_unc and cor.
    
    output:
    msg : message stating path to file.
    """
    import json 
    from numpy import array, shape
    import datetime
    
    filename = general_info['output_file'][:-3]+'json'
    
    # defining attributes  ==================================================
    if general_info['reaction'][0]=='n':
        particle = "neutron"
    elif general_info['reaction'][0]=='0':
        particle = "decay"
    else: 
        msg = "Error: projectile particle not defined for json structure."
        raise Exception(msg)   

    attributes = {'type': "DifferentialExperiment",\
                  'particle': particle, \
                  'reaction': general_info['quantity'],\
                  'nuclide': general_info['isotope'], \
                  'reaction':general_info['quantity'],\
                  'library': 'ARIADNE',\
                  'documentation':general_info['documentation'],\
                  'date':str(datetime.datetime.now().date())}
    attributes.update(features)
    # ==================================================================

    # defining data dictionary =========================================
    if general_info['quantity'] == 'nubar-prompt' or general_info['quantity'] ==  'cs' or general_info['quantity'] == 'nubar-tot':
        if shape(shape(data['einc']))[0] == 0:
            data['einc'] = array([data['einc']])
           # data['data'] = array([data['data']])
        data = {'values': (data['data']).tolist(),\
                'uncertainties': (data['rel_unc']).tolist(),\
                'correlations': (data['cor'].flatten()).tolist(),\
                'structure':{'name':'energy-in', 'type':'pointwise',\
                             'limits':(data['einc']).tolist(),\
                             'unit':data['einc_unit']},\
                'units':{'value':data['values_unit'],\
                         'uncertainty':'relative'}}
    elif general_info['quantity'] == 'PFNS':
        if shape(shape(data['einc']))[0] == 0:
            data['einc'] = array([data['einc']])
            data['eout'] = array([data['eout']])
           # data['data'] = array([data['data']])
        data = {'values': (data['data']).tolist(),\
                'uncertainties': (data['rel_unc']).tolist(),\
                'correlations': (data['cor'].flatten()).tolist(),\
                'structure':[{'name':'energy-in', 'type':'pointwise',\
                             'limits':(data['einc']).tolist(),\
                             'unit':data['einc_unit']},\
                             {'name':'energy-out', 'type':'pointwise',\
                             'limits':(data['eout']).tolist(),\
                             'unit':data['eout_unit']}],\
                'units':{'value':data['values_unit'],\
                         'uncertainty':'relative'}}
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
                  
    # ==================================================================
    
    # defining main dictionary =========================================
    main = {'attributes':attributes,'data':data}
    # ==================================================================  
       
    # write json file ===================================================
    jsonfile = open(filename,'w')
    jsonfile.write(json.dumps(main,indent=None))
    jsonfile.close()              
    # ==================================================================
    
    msg = "Output file saved in file: " + str(filename) 
    print("")       

    return msg


def write_json_output_EUCLID(general_info,data,features,identifier_iso_deriv1):
    """
    Writes json_output file for EUCLID.
    
    Input:
    general_info : dictionary containing name, isotope, quantity, reaction
    and name of the output file.
    data : dictionary containing einc_unit, eout_unit, values_unit, no_einc,
    einc, eout, values, rel_unc and cor.
    features: features related to the data.
    einc_unc: contains the identifier of the evaluated data we are calculating relative to.
    
    output:
    msg : message stating path to file.
    """
    import json 
    from numpy import array, shape, interp, zeros, arange
    import datetime
    from Manage_ReferenceData import get_reference_datacor
    from importlib import reload
    import BasicPhysicsFunctions as BPF
    reload(BPF)
    
    filename = general_info['output_file'][:-4]+'EUCLID.json'
    
    # defining attributes  ==================================================
    if general_info['reaction'][0]=='n':
        particle = "neutron"
    elif general_info['reaction'][0]=='0':
        particle = "decay"
    else: 
        msg = "Error: projectile particle not defined for json structure."
        raise Exception(msg)   


    if general_info['quantity'] == 'nubar-prompt' or general_info['quantity'] ==  'nubar-tot':
        parameter_attr = "multiplicity"
        reaction_attr = general_info['quantity']
    elif general_info['quantity'] == 'cs':
        parameter_attr = "crossSection"
        reaction_attr = general_info['reaction']
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
        
    
    attributes = {'type': "DifferentialExperiment",\
                  'particle': particle, \
                  'parameter':parameter_attr,\
                  'nuclide': general_info['isotope'], \
                  'reaction':reaction_attr,\
                  'library': 'ARIADNE',\
                  'documentation':general_info['documentation'],\
                  'date':str(datetime.datetime.now().date())}
    attributes.update(features)
    # ==================================================================

    # defining data dictionary =========================================
    expdata = {}
    if general_info['quantity'] == 'nubar-prompt' or general_info['quantity'] ==  'cs' or  general_info['quantity'] ==  'nubar-tot':
        if shape(shape(data['data']))[0] == 0:
           # expdata['einc'] = array([data['einc']])
           # expdata['data'] = array([data['data']])
            expdata = {'values': array([data['data']]).tolist(),\
                'uncertainties': data['rel_unc'].tolist(),\
                'correlations': (data['cor'].flatten()).tolist(),\
                'structure':{'name':'energy-in', 'type':'pointwise',\
                             'limits':array([data['einc']]).tolist(),\
                             'unit':data['einc_unit']},\
                'units':{'value':data['values_unit'],\
                         'uncertainty':'relative'}}
        else:
            expdata = {'values': data['data'].tolist(),\
                'uncertainties': data['rel_unc'].tolist(),\
                'correlations': (data['cor'].flatten()).tolist(),\
                'structure':{'name':'energy-in', 'type':'pointwise',\
                             'limits':(data['einc']).tolist(),\
                             'unit':data['einc_unit']},\
                'units':{'value':data['values_unit'],\
                         'uncertainty':'relative'}}
#    elif general_info['quantity'] == 'PFNS':
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
                  
    # ==================================================================

    # defining evaluated data ==========================================
    if general_info['quantity'] == 'nubar-prompt' or general_info['quantity'] == 'nubar-tot' or general_info['quantity'] == 'cs':
        # reads in nuclear data and calculates numerical derivative -------------                                                                                                                     
        # reads data ............................................................
        nucdata_lib = {'isotope':general_info['isotope'],'quantity':general_info['quantity'], \
                       'reaction':general_info['reaction'],'identifier':identifier_iso_deriv1}
        dataunccor = get_reference_datacor(nucdata_lib)
        # .......................................................................                                                                                                                     

        # converts nuclear data to units of isoptope data .......................                                                                                                                     
        if dataunccor['lattice_unit'] != data['einc_unit']:
            dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],data['einc_unit'],dataunccor['lattice_unit'])
        if dataunccor['data_unit'] != data['values_unit'] and data['values_unit'] != 'none':
            dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],data['values_unit'],dataunccor['data_unit'])
        # .......................................................................                                           

        # interpolates nuclear data to lattice of experiment                                                                                                                                                    
        values = array(interp(data['einc'],dataunccor['lattice'],dataunccor['data']))
        # -----------------------------------------------------------------------
        evaluateddata = {'library':identifier_iso_deriv1,'values':values.tolist()}
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
    # ==================================================================

    # defining sensitivities ===========================================
    energy = [1e-11, 2.9694e-09, 1.0364e-08, 2.999e-08, 4.0991e-08, 4.9445e-08,\
          7.1941e-08, 1.0467e-07, 1.52e-07, 2.2159e-07, 2.511e-07, 2.8453e-07,\
          3.2242e-07, 3.6535e-07, 4.14e-07, 6.4121e-07, 1.13e-06, 3.06e-06,\
          8.32e-06, 2.26e-05, 6.14e-05, 0.000167, 0.000454, 0.001235, 0.00335,\
          0.00912, 0.0248, 0.0676, 0.184, 0.303, 0.5, 0.823, 1.353, 1.738, 2.232,\
          2.865, 3.68, 4.7237, 5.3526, 6.07, 6.8729, 7.79, 8.825, 10.0, 12.0, 13.0,\
          13.5, 14.0, 14.5, 15.0, 17.0, 20.0]
    energyunit = "MeV"

    # check if same units ..............................................
    if energyunit != data['einc_unit']:
        msg = "Wrong unit of data."
        raise Exception(msg)

    if shape(shape(data['einc']))[0] > 0:
        sens = zeros([shape(data['einc'])[0],shape(energy)[0]],dtype=float)

        for index1 in arange(0,shape(data['einc'])[0]):
            eincExp = data['einc'][index1]
            for index2 in arange(0,shape(energy)[0]-1):
                if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                    sens[index1,index2] = 1.0
    else:
        sens = zeros(shape(energy)[0],dtype=float)
        eincExp = data['einc']
        for index2 in arange(0,shape(energy)[0]-1):
            if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                sens[index2] = 1.0

                     
    sensitivities = {'values':sens.tolist(),\
                     'structure':{'name':'energy-column','type':'histogram',\
                                  'limits':energy,'unit':energyunit}}
    # ==================================================================
    
    # defining main dictionary =========================================
    main = {'attributes':attributes,'experimental data':expdata,\
            'evaluated data':evaluateddata,'sensitivities':sensitivities}
    # ==================================================================  
       
    # write json file ===================================================
    jsonfile = open(filename,'w')
    jsonfile.write(json.dumps(main,indent=None))
    jsonfile.close()              
    # ==================================================================
    
    msg = "Output file saved in file: " + str(filename) 
    print("")       

    return msg

def write_json_output_ratio_EUCLID(general_info,data,features,identifier_iso_deriv1,reference):
    """
    Writes json_output file for EUCLID.
    
    Input:
    general_info : dictionary containing name, isotope, quantity, reaction
    and name of the output file.
    data : dictionary containing einc_unit, eout_unit, values_unit, no_einc,
    einc, eout, values, rel_unc and cor.
    features: features related to the data.
    einc_unc: contains the identifier of the evaluated data we are calculating relative to.
    
    output:
    msg : message stating path to file.
    """
    import json 
    from numpy import array, shape, interp, zeros, arange
    import datetime
    from Manage_ReferenceData import get_reference_datacor
    from importlib import reload
    import BasicPhysicsFunctions as BPF
    reload(BPF)
    
    filename = general_info['output_file'][:-4]+'ratiodataEUCLID.json'
    
    # defining attributes  ==================================================
    if general_info['reaction'][0]=='n':
        particle = "neutron"
    elif general_info['reaction'][0]=='0':
        particle = "decay"
    else: 
        msg = "Error: projectile particle not defined for json structure."
        raise Exception(msg)   

    if general_info['quantity'] == 'cs' and reference['quantity'] == 'cs':
        parameter_attr = "crossSection/crossSection"
        reaction_attr = general_info['reaction']+"/"+reference['reaction']
    else:
        msg = "Combination of ratio quantities not yet defined for json structure."
        raise Exception(msg) 

    attributes = {'type': "DifferentialExperiment",\
                  'particle': particle, \
                  'parameter':parameter_attr,\
                  'nuclide': general_info['isotope']+"/"+reference['isotope'], \
                  'reaction':reaction_attr,\
                  'library': 'ARIADNE',\
                  'documentation':general_info['documentation'],\
                  'date':str(datetime.datetime.now().date())}
    attributes.update(features)
    # ==================================================================

    # defining data dictionary =========================================
    expdata = {}
    if general_info['quantity'] == 'nubar-prompt' or  general_info['quantity'] ==  'cs' or general_info['quantity'] ==  'nubar-tot':
        if shape(shape(data['data']))[0] == 0:
            expdata['einc'] = array([data['einc']])
        expdata = {'values': data['ratio'].tolist(),\
                'uncertainties': data['relunc_ratio'].tolist(),\
                'correlations': (data['cor_ratio'].flatten()).tolist(),\
                'structure':{'name':'energy-in', 'type':'pointwise',\
                             'limits':(data['einc']).tolist(),\
                             'unit':data['einc_unit']},\
                'units':{'value':data['ratio_unit'],\
                         'uncertainty':'relative'}}
#    elif general_info['quantity'] == 'PFNS':
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
                  
    # ==================================================================

    # defining evaluated data ==========================================
    if general_info['quantity'] == 'nubar-prompt' or general_info['quantity'] == 'nubar-tot' or general_info['quantity'] == 'cs':
        # reads in nuclear data and calculates numerical derivative -------------                                                                                                                     
        # reads data ............................................................
        # main reaction
        nucdata_lib = {'isotope':general_info['isotope'],'quantity':general_info['quantity'], \
                       'reaction':general_info['reaction'],'identifier':identifier_iso_deriv1}
        dataunccor = get_reference_datacor(nucdata_lib)

        # reference
        nucdata_lib_ref =  {'isotope':reference['isotope'],'quantity':reference['quantity'], \
                           'reaction':reference['reaction'],'identifier':reference['identifier']}
        dataunccor_ref = get_reference_datacor(nucdata_lib_ref)
        # .......................................................................                                                                                                                     

        # converts nuclear data to units of isoptope data .......................  
        # main reaction                                                                                                                   
        if dataunccor['lattice_unit'] != data['einc_unit']:
            dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],data['einc_unit'],dataunccor['lattice_unit'])

        if dataunccor['data_unit'] != data['values_unit']:
            dataunccor['data'] = BPF.conversion_to_newUnits(dataunccor['data'],data['values_unit'],dataunccor['data_unit'])

        # reference                                                                                                             
        if dataunccor_ref['lattice_unit'] != data['einc_unit']:
            dataunccor_ref['lattice'] = BPF.conversion_to_newUnits(dataunccor_ref['lattice'],data['einc_unit'],dataunccor_ref['lattice_unit'])

        if dataunccor_ref['data_unit'] != data['values_unit']:
            dataunccor_ref['data'] = BPF.conversion_to_newUnits(dataunccor_ref['data'],data['values_unit'],dataunccor_ref['data_unit'])
        # .......................................................................                                           

        # interpolates nuclear data to lattice of experiment         
        # main reaction
        values = array(interp(data['einc'],dataunccor['lattice'],dataunccor['data']))
        
        # reference   
        values_ref = array(interp(data['einc'],dataunccor_ref['lattice'],dataunccor_ref['data']))

        ratios_eval = values/values_ref
        # -----------------------------------------------------------------------
        evaluateddata = {'library':identifier_iso_deriv1,'values':ratios_eval.tolist()}
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg)   
    # ==================================================================

    # defining sensitivities ===========================================
    energy = [1e-11, 2.9694e-09, 1.0364e-08, 2.999e-08, 4.0991e-08, 4.9445e-08,\
          7.1941e-08, 1.0467e-07, 1.52e-07, 2.2159e-07, 2.511e-07, 2.8453e-07,\
          3.2242e-07, 3.6535e-07, 4.14e-07, 6.4121e-07, 1.13e-06, 3.06e-06,\
          8.32e-06, 2.26e-05, 6.14e-05, 0.000167, 0.000454, 0.001235, 0.00335,\
          0.00912, 0.0248, 0.0676, 0.184, 0.303, 0.5, 0.823, 1.353, 1.738, 2.232,\
          2.865, 3.68, 4.7237, 5.3526, 6.07, 6.8729, 7.79, 8.825, 10.0, 12.0, 13.0,\
          13.5, 14.0, 14.5, 15.0, 17.0, 20.0]
    energyunit = "MeV"

    # check if same units ..............................................
    if energyunit != data['einc_unit']:
        msg = "Wrong unit of data."
        raise Exception(msg)

    # main reaction
    if shape(shape(data['einc']))[0] > 0:
        sens = zeros([shape(data['einc'])[0],shape(energy)[0]],dtype=float)

        for index1 in arange(0,shape(data['einc'])[0]):
            eincExp = data['einc'][index1]
            for index2 in arange(0,shape(energy)[0]-1):
                if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                    sens[index1,index2] = 1.0/values_ref[index1]
    else:
        sens = zeros(shape(energy)[0],dtype=float)
        eincExp = data['einc']
        for index2 in arange(0,shape(energy)[0]-1):
            if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                sens[index2] = 1.0/values_ref
                
    if  general_info['quantity'] == 'cs':
        parameter_attr = "crossSection"
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg) 

    attributes_sens_main = {'parameter':parameter_attr,\
                  'nuclide': general_info['isotope'],\
                  'reaction':general_info['reaction']}
    
    sensitivities_main = {'values':sens.tolist(),\
                     'structure':{'name':'energy-column','type':'histogram',\
                                  'limits':energy,'unit':energyunit},\
                          'attributes':attributes_sens_main}

    # reference reaction
    if shape(shape(data['einc']))[0] > 0:
        sens = zeros([shape(data['einc'])[0],shape(energy)[0]],dtype=float)

        for index1 in arange(0,shape(data['einc'])[0]):
            eincExp = data['einc'][index1]
            for index2 in arange(0,shape(energy)[0]-1):
                if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                    sens[index1,index2] = -1.0*values[index1]/(values_ref[index1]*values_ref[index1])
    else:
        sens = zeros(shape(energy)[0],dtype=float)
        eincExp = data['einc']
        for index2 in arange(0,shape(energy)[0]-1):
            if energy[index2] <= eincExp and energy[index2+1] > eincExp:
                sens[index2] = -1.0*values/(values_ref*values_ref)

    if  reference['quantity'] == 'cs':
        parameter_attr = "crossSection"
    else:
        msg = "Quantity not yet defined for json structure."
        raise Exception(msg) 

    attributes_sens_ref = {'parameter':parameter_attr,\
                  'nuclide': reference['isotope'],\
                  'reaction':reference['reaction']}
    
    sensitivities_ref = {'values':sens.tolist(),\
                     'structure':{'name':'energy-column','type':'histogram',\
                                  'limits':energy,'unit':energyunit},\
                         'attributes':attributes_sens_ref}
                

    sensitivities = {'sensitivity main observable':sensitivities_main,\
                     'sensitivity reference observable':sensitivities_ref}
    # ==================================================================
    
    # defining main dictionary =========================================
    main = {'attributes':attributes,'experimental data':expdata,\
            'evaluated data':evaluateddata,'sensitivities':sensitivities}
    # ==================================================================  
       
    # write json file ===================================================
    jsonfile = open(filename,'w')
    jsonfile.write(json.dumps(main,indent=None))
    jsonfile.close()              
    # ==================================================================
    
    msg = "Output file saved in file: " + str(filename) 
    print("")       

    return msg

def test_write_xml_output(): 
    """ 
    very simple test for test_write_xml_output()
    """
    from numpy import array
    

    general_info = {'name': 'Vorobyev', 'isotope': 'U-235', 'quantity': 'PFNS', 'reaction': 'n,f', \
                'documentation' : "this is only a test.",\
                'output_file' : BASE_DIR+'Example/PFNS_test.xml'}
    data_test = {'einc_unit':'MeV','eout_unit': 'MeV','values_unit': "none",\
             'no_einc': 2,'einc': [0.5,0.5,1.0,1.0],'eout': array([1,2,1,2]),'data':[10,20,0.1,0.2], \
             'rel_unc': [1.5,1.25,10.0,12.0] ,\
             'cor' : array([[1,0.2,0.3,0.4],[0.2, 1, 0.5, 0.6], \
                      [0.3,0.5,1.0,0.7],[0.4,0.6,0.7,1]])}

    write_xml_output(general_info,data_test)

    return None
