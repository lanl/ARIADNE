# This module contains a GLS evaluation code for PFNS.

from numpy import isscalar, array, shape, zeros, arange, loadtxt, ndarray, matrix, linalg, dot, inner, sqrt, diag, interp


def Evaluate(ExpInput,PriorInput,homedir,CommonGrid=False,PPPCorrection=False,Plot=False,MaxwRatioEval=False,LogEval=False,OutputFormat="Dece"):
    """
    Provides evaluated data and covariances using the GLS algorithm.
    
    
    Input parameters:
    ---------------------------------------------
    ExpInput: dictonary of experimental inputs containing
              flag whether this is a GMA or ARIADNE file.
              The second argument is a dictionary pointing 
              to the experimental information
    PriorInput: dictionary with one array on observables.
              One array with file descriptors of mv, relunc and cov.
    CommonGrid: default False. If True, exp. data will be put on 
              common grid of Prior before evaluation.
    PPPCorrection: default False. If True, the Chiba-Smith 
              PPP correction will be enforced on the evaluation.
    Plot: default False, If True, evluated results will
              be plotted compared to experimental data and prior.
    OuputFormat: default Dece output format.
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    Evaluated data in output files and plots if selected.
    ---------------------------------------------
    """
    from numpy import sum, concatenate

    # reading input -----------------------------
    if ExpInput['Code'] == "ARIADNE":
        ExpData = ReadAriadneInput(ExpInput['ReadOption'],ExpInput['Data'])
    elif ExpInput['Code'] == "GMA":
        ExpData = ReadGMAInput(ExpInput['data'])
    
    PriorData = ReadPriorInput(PriorInput)

    ExpData = CutExpToPriorGird(ExpData,PriorData)
    # ---------------------------------------------

    # evaluation process --------------------------
    if CommonGrid:
        ExpData = BringToCommonGrid(ExpData,PriorData)
    
    Sens = CalculateSensitivities(ExpData,PriorData)

    EvaluatedData = GLS(ExpData,PriorData,Sens,PPPCorrection,MaxwRatioEval,LogEval)
    # ---------------------------------------------

    # output formatting ---------------------------
    if Plot:
        Plot(ExpData,PriorData,EvaluatedData)

    OutputFilename = WriteEvaluationOutput(OutputFormat,homedir,PriorInput,EvaluatedData)

    print("The evaluation output is saved in ", OutputFilename)
    # ---------------------------------------------
    print("Evaluate executed successfully.")
    return OutputFilename


###############################################################
def ReadAriadneInput(ReadOption,ExpInput):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpInput: A dictionary pointing to  the experimental information
              is given.
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    ExpData: 
    ---------------------------------------------
    """
    import json
    import itertools
    from numpy import sum, concatenate
    
    # initialization ------------------------------------------
    Einc = []
    Mv = []
    RelUnc = []
    observable = []
    # ---------------------------------------------------------

    # reading mv, Einc, observable  ---------------------------
    if ReadOption == 'Ariadneclass':
        for Exp in ExpInput:
            Einc.append(Exp.einc)
            Mv.append(Exp.data_total['data'])
            RelUnc.append(Exp.data_total['rel_unc'])
            
            for index2 in arange(0,Exp.no_val):
                observable.append(Exp.isotope+Exp.general_info['reaction']+Exp.general_info['quantity'])
    elif ReadOption == 'json':
        # opening file, getting name, einc and no of elements -----------------
        with open (ExpInput) as jsonfile:
            experimentdoc = json.load( jsonfile )

        print('------------------------------------')    
        print('Experiments included in evaluation:')
        for experiment in experimentdoc:
            print(experiment['attributes']['Author'])
            Einc   = concatenate([Einc, experiment['data']['structure'][1]['limits'][0]])
            Mv     = concatenate([Mv, experiment['data']['values']])
            RelUnc = concatenate([RelUnc,experiment['data']['uncertainties']])
            
            for index2 in arange(0,shape(experiment['data']['structure'][1]['limits'])[0]):   
                observable.append(experiment['attributes']['nuclide']\
                                  +experiment['attributes']['reaction'])
            
        print('------------------------------------') 
    # ---------------------------------------------------------------------
   
        
    else:
        print("Warning: no experimental data read. ReadOption ", ReadOption, " not implemented.")

        
    Einc_arr = array(Einc)
    Mv_arr = array(Mv)
    RelUnc_arr = array(RelUnc)
    observables_arr = array(observable)
    # ---------------------------------------------------------
    
    # reading cor ---------------------------------------------
    no_data = max(shape(Einc_arr))
    cor = zeros([no_data,no_data],dtype=float)
    cov = zeros([no_data,no_data],dtype=float)
    startindex = 0
    
    if ReadOption == 'Ariadneclass':
        for Exp in ExpInput:
            endindex = startindex+Exp.no_val
            cor[startindex:endindex,startindex:endindex]=Exp.data_total['cor']
            startindex = endindex
    elif ReadOption == 'json':
        with open (ExpInput) as jsonfile:
            experimentdoc = json.load( jsonfile )

        for experiment in experimentdoc:
            corflat =  experiment['data']['correlations']
            endindex = startindex+shape(experiment['data']['structure'][1]['limits'][0])[0]
                
            for index1 in arange(startindex,endindex):
                for index2 in arange(startindex,endindex):
                    cor[index1,index2] = corflat[(index1-startindex)*(endindex-startindex)+index2-startindex]
                    
            startindex = endindex
        
    else:
        print("Warning: no experimental data read. ReadOption ", ReadOption, " not implemented.")        
    # ---------------------------------------------------------
    import matplotlib.pyplot as plt
    plt.pcolormesh(arange(0,shape(cor)[0]),arange(0,shape(cor)[0]),cor,vmin=-1, vmax=1,cmap='bwr')
    plt.xlabel("Number of data points")
    plt.ylabel("Number of data points")
    plt.colorbar()
    
    plt.figure(2)
    plt.xscale('log')
    plt.xlabel("Outgoing Neutron Energy (MeV)")
    plt.ylabel("PFNS (1/ MeV)")
    plt.errorbar(Einc_arr,Mv_arr,RelUnc_arr*Mv_arr,fmt='o')

    # preparing output ----------------------------------------
    for index1 in arange(0,no_data):
        for index2 in arange(0,no_data):
            cov[index1,index2] = cor[index1,index2]*RelUnc_arr[index1]*RelUnc_arr[index2]*Mv_arr[index1]*Mv_arr[index2]

    ExpData = {'observables':observables_arr,'Einc':Einc_arr,'Mv':Mv_arr,'cov':cov}
    # ---------------------------------------------------------
    return ExpData


###############################################################

def CutExpToPriorGird(ExpData,PriorData):
    """
    Cuts experimental data data prior grid.

    
    Input parameters:
    ---------------------------------------------
    ExpData: A dictionary pointing to  the experimental information
              is given.
    PriorData: Prior data.
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    ExpData: A dictionary pointing to  the experimental information
              is given. In this dictionary, the experimental data are
              cut to be within the energy grid of the prior.
    ---------------------------------------------
    """

    # Initialization ----------------------------
    PriorEinc = PriorData['Einc']
    ExpEinc = ExpData['Einc']
    ExpMv = ExpData['Mv']
    ExpCov = ExpData['cov']

    maxEinc = max(PriorEinc)
    count = 0
    for E in ExpEinc:
        if E > maxEinc:
            count = count+1

            
    if count > 0:
        dim = shape(ExpEinc)[0]-count
        ExpEincNew = zeros(dim,dtype=float)
        ExpMvNew = zeros(dim,dtype=float)
        ExpCovNew = zeros([dim,dim],dtype=float)
    # ---------------------------------------------

    
    # Cutting data --------------------------------    
        count = 0
        for index in arange(0,shape(ExpEinc)[0]):
            if ExpEinc[index] <= maxEinc:
                ExpEincNew[count]= ExpEinc[index]
                ExpMvNew[count] = ExpMv[index]
                count = count+1

        count1 = 0
        count2 = 0
        for index1 in arange(0,shape(ExpEinc)[0]):
            if ExpEinc[index1] <= maxEinc:
                for index2 in arange(0,shape(ExpEinc)[0]):
                    if ExpEinc[index2] <= maxEinc:
                        ExpCovNew[count1,count2] = ExpCov[index1,index2]
                        count2 = count2+1
                count2 = 0
                count1 = count1 + 1
    # --------------------------------------------


    # Preparing output ---------------------------    
        ExpData['Einc'] = ExpEincNew
        ExpData['Mv'] = ExpMvNew
        ExpData['cov'] = ExpCovNew
    # --------------------------------------------
        print("\n Some experimental data not within prior grid.")
        print("-------------------------------------------")
    else:
        print("\n All experimental data within prior grid.")
        print("-------------------------------------------")
	
    return ExpData
	
###############################################################


def ReadGMAInput(ExpInput):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpInput: A dictionary pointing to  the experimental information
              is given.
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    ExpData: 
    ---------------------------------------------
    """
    ExpData = 0
    print("ReadGMAInput: Not yet implemented")

    
    return ExpData


###############################################################


def ReadPriorInput(PriorInput):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    PriorInput: A dictionary pointing to  the experimental information
              is given.
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    PriorData: 
    ---------------------------------------------
    """
    
    # initialization ------------------------------------------
    PriorData = 0

    Einc = []
    Mv = []
    RelUnc = []
    observable = []
    
    general_info = PriorInput['general_info']
    mv_files = PriorInput['Data']['mv']
    cor_files = PriorInput['Data']['cor']
    no_obs = shape(general_info['isotope'])[0]
    # ---------------------------------------------------------
    
    # reading mv, Einc, observable  ---------------------------
    for index in arange(0,no_obs):
        data = loadtxt(mv_files[index],comments='#')
        
        Einc.append(data[:,0])
        Mv.append(data[:,1])
        RelUnc.append(data[:,2])

        for index2 in arange(0,shape(data)[0]):
            observable.append(general_info['isotope'][index]+general_info['reaction'][index]+general_info['quantity'][index])

    Einc_arr = array(Einc).flatten()
    Mv_arr = array(Mv).flatten()
    RelUnc_arr = array(RelUnc).flatten()
    observables_arr = array(observable)
    # ---------------------------------------------------------
    
    # reading cor ---------------------------------------------
    no_data = max(shape(observables_arr))
    cor = zeros([no_data,no_data],dtype=float)
    cov = zeros([no_data,no_data],dtype=float)
    startindex = 0

    for index in arange(0,no_obs):
        data = loadtxt(cor_files[index],comments='#')
        endindex = startindex+shape(data)[0]
        cor[startindex:endindex,startindex:endindex]=data
        startindex = endindex
    # ---------------------------------------------------------

    # preparing output ----------------------------------------
    for index1 in arange(0,no_data):
        for index2 in arange(0,no_data):
            cov[index1,index2] = cor[index1,index2]*RelUnc_arr[index1]*RelUnc_arr[index2]*Mv_arr[index1]*Mv_arr[index2]/(100.0**2.0)

    PriorData = {'observables':observables_arr,'Einc':Einc_arr,'Mv':Mv_arr,'cov':cov}
    # ---------------------------------------------------------

    return PriorData


###############################################################


def CommonGrid(ExpData,PriorData):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpData: either from ReadAriadneInput or ReadGMAInput
            contains experimental data, cov, observable and shape information
    PriorData: from ReadPriorInput
            contains prior data, cov and observable information
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    ExpDataCommonGrid: now on common grid of prior for a specific observable. 
    ---------------------------------------------
    """
    ExpDataCommonGrid = 0
    print("CommonGrid: Not yet implemented")
    
    return ExpDataCommonGrid


###############################################################



def CalculateSensitivities(ExpData,PriorData):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpData: either from ReadAriadneInput or ReadGMAInput
            contains experimental data, cov, observable and shape information
    PriorData: from ReadPriorInput
            contains prior data, cov and observable information
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    Sens: Sensitivity profiles
    ---------------------------------------------
    """
    
    # initialization ------------------------------------------
    Einc_prior = PriorData['Einc']
    Einc_exp = ExpData['Einc']
    Sens = zeros([shape(Einc_exp)[0],shape(Einc_prior)[0]],dtype=float)
    # ---------------------------------------------------------
    
    # initialization ------------------------------------------
    for (index1,Ee) in enumerate(Einc_exp):
        for (index2,Ep) in enumerate(Einc_prior):
  #          if ExpData['observables'][index1] == PriorData['observables'][index2]:
            if (Ee >= Einc_prior[index2-1]) and (Ee < Ep) and index2>0:
                Sens[index1,index2] = (Ee-Einc_prior[index2-1])/(Ep-Einc_prior[index2-1])
            if (Ee >= Ep) and (Ee < Einc_prior[index2+1]) and index2<shape(Einc_prior)[0]:
                Sens[index1,index2] = (Ee-Einc_prior[index2+1])/(Ep-Einc_prior[index2+1])
            #else:
            #    "Warning: different observables not yet implemented. Sens set to 0."
    # ---------------------------------------------------------
    return Sens


###############################################################


def GLS(ExpData,PriorData,Sens,PPPCorrection,MaxwRatioEval,LogEval):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpData: either from ReadAriadneInput or ReadGMAInput
            contains experimental data, cov, observable and shape information
    PriorData: from ReadPriorInput
            contains prior data, cov and observable information
    Sens: Sensitivity profiles
    PPPCorrection: 
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    EvaluatedData
    ---------------------------------------------
    """
    from numpy import log, exp
    
    PriorEinc = PriorData['Einc']
    PriorMv = PriorData['Mv']
    PriorCov = PriorData['cov']
    ExpEinc = ExpData['Einc']
    ExpMv = ExpData['Mv']
    ExpCov = ExpData['cov']


    # Conversion to ratio to Maxwellian ---------------------------------------------------------------------------
    if MaxwRatioEval:
        T = 1.42 # MeV
        maxwPri = maxwt(PriorEinc,T)
        PriorMv = PriorMv/maxwPri
        maxwExp = maxwt(ExpEinc,T)
        ExpMv   = ExpMv/maxwExp

        for index1 in arange(0,shape(ExpMv)[0]):
            for index2 in arange(0,shape(ExpMv)[0]):
                ExpCov[index1,index2] = ExpCov[index1,index2]/(maxwExp[index1]*maxwExp[index2])
                
        for index1 in arange(0,shape(PriorMv)[0]):
            for index2 in arange(0,shape(PriorMv)[0]):
                PriorCov[index1,index2] = PriorCov[index1,index2]/(maxwPri[index1]*maxwPri[index2])
    # -----------------------------------------------------------------------------------------------------------

    # Conversion to log space -----------------------------------------------------------------------------------
    if LogEval:
        for index1 in arange(0,shape(ExpMv)[0]):
            for index2 in arange(0,shape(ExpMv)[0]):
                ExpCov[index1,index2] = ExpCov[index1,index2]/(ExpMv[index1]*ExpMv[index2])
        ExpMv = log(ExpMv)

        for index1 in arange(0,shape(PriorMv)[0]):
            for index2 in arange(0,shape(PriorMv)[0]):
                PriorCov[index1,index2] = PriorCov[index1,index2]/(PriorMv[index1]*PriorMv[index2])
        PriorMv = log(PriorMv)
    # -----------------------------------------------------------------------------------------------------------
    
    if PPPCorrection:
        prior_int = array(interp(ExpEinc,PriorEinc,PriorMv))
        for index1 in arange(0,shape(ExpMv)[0]):
            for index2 in arange(0,shape(ExpMv)[0]):
                ExpCov[index1,index2] = ExpCov[index1,index2]*prior_int[index1]*prior_int[index2]/(ExpMv[index1]*ExpMv[index2])
        iterations = 10

    # Evaluation ---------------------------------
    EvalMv = zeros(shape(PriorMv)[0],dtype=float)
    EvalCov = zeros([shape(PriorMv)[0],shape(PriorMv)[0]],dtype=float)

    Qinv = linalg.inv(dot(dot(Sens,PriorCov),Sens.transpose())+ExpCov)

    EvalCov = PriorCov - dot(dot(dot(dot(PriorCov,Sens.transpose()),Qinv),Sens),PriorCov)
    d = ExpMv-dot(Sens,PriorMv)
    EvalMv = PriorMv+dot(dot(dot(EvalCov,Sens.transpose()),linalg.inv(ExpCov)),d)   
    
    # PPP iterations ................................    
    if PPPCorrection:
        # executing Chiba-Smith correction.
        for iter in range(1,iterations):
            prior_intNew = array(interp(ExpEinc,PriorEinc,EvalMv))
            print('PPP execution, deviation of previous iteration from current iteration prior (should converge with number of iterations):', prior_intNew/prior_int)
        
            for index1 in arange(0,shape(ExpMv)[0]):
                for index2 in arange(0,shape(ExpMv)[0]):
                    ExpCov[index1,index2] = ExpCov[index1,index2]*prior_intNew[index1]*prior_intNew[index2]/(prior_int[index1]*prior_int[index2])

            Qinv = linalg.inv(dot(dot(Sens,PriorCov),Sens.transpose())+ExpCov)
            
            EvalCovPar = PriorCov - dot(dot(dot(dot(PriorCov,Sens.transpose()),Qinv),Sens),PriorCov)
            d = ExpMv-dot(Sens,PriorMv)
            EvalMv = PriorMv+dot(dot(dot(EvalCov,Sens.transpose()),linalg.inv(ExpCov)),d)            
            
            prior_int = prior_intNew
    # ---------------------------------------------


    # Chi^2 ---------------------------------------
    noofexp = shape(ExpMv)[0] 
    chi2    = dot(dot(d.transpose(),Qinv),d)/noofexp
    print('Chi^2:', chi2)
    # ---------------------------------------------


    # Conversion from ratio to Maxwellian ---------------------------------------------------------------------------
    if MaxwRatioEval:
        print("Tranforming back from ratio to Maxwellian")
        maxwEv   = maxwt(PriorEinc,T)
        EvalMv   = EvalMv*maxwEv
        
                
        for index1 in arange(0,shape(PriorMv)[0]):
            for index2 in arange(0,shape(PriorMv)[0]):
                EvalCov[index1,index2] = EvalCov[index1,index2]*(maxwEv[index1]*maxwEv[index2])
    # -----------------------------------------------------------------------------------------------------------

    # Conversion to log space -----------------------------------------------------------------------------------
    if LogEval:
        EvalMv = exp(EvalMv)
        for index1 in arange(0,shape(EvalMv)[0]):
            for index2 in arange(0,shape(EvalMv)[0]):
                EvalCov[index1,index2] = EvalCov[index1,index2]*(EvalMv[index1]*EvalMv[index2])
    # -----------------------------------------------------------------------------------------------------------

    
    
    EvalData = {'observables':PriorData['observables'],'Einc':PriorData['Einc'],'Mv':EvalMv,'cov':EvalCov}
    
    return EvalData


###############################################################

def maxwt(E,T):
    from numpy import sqrt, exp,zeros, max, shape
    pi = 3.14159265359
    
    maxwt = zeros(shape(E)[0],dtype=float)
    for index in arange(0,shape(E)[0]):
        maxwt[index] = sqrt(E[index])*exp(-E[index]/T)*(2./sqrt(pi)/sqrt(T)/T)
                     
    return maxwt
    

def Plot(ExpData,PriorData,EvaluatedData):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    ExpData: either from ReadAriadneInput or ReadGMAInput
            contains experimental data, cov, observable and shape information
    PriorData: from ReadPriorInput
            contains prior data, cov and observable information
    EvaluatedData
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    EvaluatedData
    ---------------------------------------------
    """
    print("Plot: Not yet implemented")
    
    return 0


###############################################################

def WriteEvaluationOutput(OutputFormat,homedir,PriorInput,EvalData):
    """
    Provides 
    
    
    Input parameters:
    ---------------------------------------------
    EvaluatedData
    ---------------------------------------------
    
    Returns:
    ---------------------------------------------
    OutputFilename
    ---------------------------------------------
    """
    
    # initialization ----------------------------
    filenamemv = 'GLS.dat'
    filenamecor = 'GLS.dat'
    general_info = PriorInput['general_info']
    no_obs = shape(general_info['isotope'])[0]
    index0 = 0
    dim = shape(array(EvalData['Einc']))[0]
    cor = zeros(shape(array(EvalData['cov'])),dtype=float)
    RelUnc = sqrt(diag(EvalData['cov']))
    # ---------------------------------------------
    
    # getting correlation matrix ----------------
    for index1 in arange(shape(EvalData['cov'])[0]):
        for index2 in arange(shape(EvalData['cov'])[0]):
            cor[index1,index2] = EvalData['cov'][index1,index2]/(RelUnc[index1]*RelUnc[index2])
    # ---------------------------------------------
    
    # writing to mv, cor file --------------------
    for indexobs in arange(no_obs):
       observable = general_info['isotope'][indexobs]+general_info['reaction'][indexobs]+general_info['quantity'][indexobs]

    # ARIADNE txt file format --------------------
       if OutputFormat == 'ARIADNEtxt':
           
        # ARIADNE mv header ----------------------
           filenamemv = homedir + 'Eval_MvUnc'+observable+'GLS.dat'
           fmv = open(filenamemv,'w')
           string='#'+observable+'\n'
           fmv.write(string)
           string='# Einc( MeV)   mv (b) RelUnc(%)\n'
           fmv.write(string)
        # ----------------------------------------
        
        # ARIADNE cor header ---------------------
           filenamecor = homedir + 'Eval_Cor'+observable+'GLS.dat'
           fcor = open(filenamecor,'w')
           string='#'+observable+' cor values\n'
           fcor.write(string)
        # ----------------------------------------
    
        # ARIADNE mv file ------------------------
           for (index1,observableE) in enumerate(EvalData['observables'][index0:]):
               if observableE == observable:
                   string = str(EvalData['Einc'][index1])+'\t'+str(EvalData['Mv'][index1])+'\t'+\
                       str(RelUnc[index1]*100.0/EvalData['Mv'][index1])+'\n'
                   fmv.write(string)
               else:
                   break
        # ----------------------------------------
        
        # ARIADNE cor file  ----------------------
           for indexcor in arange(index0,index1+1):
               string=' '.join(map(str,cor[indexcor,index0:index1+1]))+'\n'
               fcor.write(string)
               
        # ----------------------------------------
        
        # closing --------------------------------
           index0 = index1
           fmv.close()
           fcor.close()    
        # ----------------------------------------
    # ---------------------------------------------
    # ---------------------------------------------

    return [filenamemv,filenamecor]


###############################################################
