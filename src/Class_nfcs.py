# this module contains the classes for fission cross sections experiments

class nfcs(object):
    # mother class nfcs
    
    def __init__(self,general_info,data,unc_iso,features):
        from numpy import shape, zeros
        from BasicPhysicsFunctions import conversion_to_newUnits

        # general information ---------------------------------------------------------
        self.general_info = general_info
        # WIP: isotope checking
        self.isotope = general_info['isotope']
        # WIP: error message if not n,f
        # WIP: error message if n,f and Einc = 0
        self.features = features
        # ----------------------------------------------------------------------------
        
        # loading cross section values and their incident neutron energies -----------
        self.values = data['values']
        self.values_unit = data['values_unit']
        if shape(shape(self.values))[0] == 0:
            self.no_val = 1
        else:
            self.no_val = shape(self.values)[0]
        
        # comment: number of Einc is not needed for cross sections
        self.einc = data['einc']
        self.einc_unit = data['einc_unit']
        # ---------------------------------------------------------------------------
        
        # energy resolution ---------------------------------------------------------

        self.edges_for_chw_integrated_derivative = None
        self.tof_length_in_resolution = False
        self.calculate_energy_resolution = False

        if 'resolution' in unc_iso['einc_unc']:
            resolution = unc_iso['einc_unc']['resolution']

            if 'tof_length' in resolution:
                self.tof_length_val = resolution['tof_length']
                self.tof_length_unit = resolution['tof_length_unit']
                self.tof_length_in_resolution = True

            if 'chw_edges' in resolution:
                assert 'chw_edges_unit' in resolution, ("chw_edges given but no chw_edges_unit")

                if resolution['chw_edges_unit'] != self.einc_unit:
                    resolution['chw_edges'] = conversion_to_newUnits(resolution['chw_edges'], self.einc_unit, resolution['chw_edges_unit'])

                self.edges_for_chw_integrated_derivative = resolution['chw_edges']
            
            if 'energy_resolution' in resolution:
                self.calculate_energy_resolution = True

                if resolution['energy_resolution_unit'] != self.einc_unit:
                    resolution['energy_resolution'] = conversion_to_newUnits(resolution['energy_resolution'], self.einc_unit, resolution['energy_resolution_unit'])
            
                self.energy_resolution = resolution['energy_resolution']

                # warning for low energy resolution broadening 
                einc_eV = conversion_to_newUnits(self.einc, 'eV', self.einc_unit)
                if len(shape(einc_eV)) < 1:
                    if einc_eV<1e-2:
                        raise ValueError("Energy resolution broadening for regions below 100eV has not be validated")
                else:    
                    if min(einc_eV)<1e-2:
                        raise ValueError("Energy resolution broadening for regions below 100eV has not be validated")

        ### Check for old inputs
        if 'en_rsl' in unc_iso['einc_unc']:
            raise ValueError("en_rsl is no longer a valid arguement - calculate_cov_eunc_legacy() functions were maintained in ToFCovariance module and can be recovered if necessary")

        if 'trsl' in unc_iso['einc_unc']:
            raise ValueError("trsl is no longer a valid arguement - calculate_cov_eunc_legacy() functions were maintained in ToFCovariance module and can be recovered if necessary")

        # energy uncertainties ------------------------------------------------------
        if 'en_err' in unc_iso['einc_unc']:
            self.do_energy_unc = True
            self.enerr_unc = unc_iso['einc_unc']['en_err']['enerr_unc']
            self.enerr_unc_unit = unc_iso['einc_unc']['en_err']['enerr_unc_unit']
            self.enerr_unc_type = unc_iso['einc_unc']['en_err']['enerr_unc_type']
            self.enerr_unc_type_arg = unc_iso['einc_unc']['en_err']['enerr_unc_type_arg']
        else:
            self.do_energy_unc = False

        if 'tof_length' in unc_iso['einc_unc']:
            self.do_tof_length_unc = True

            if self.tof_length_in_resolution:
                assert(self.tof_length_val == unc_iso['einc_unc']['tof_length']['value']), "tof_length in resolution is different from that in einc_unc"
                assert(self.tof_length_unit == unc_iso['einc_unc']['tof_length']['value_unit']), "tof_length in resolution is different from that in einc_unc"
            else:
                self.tof_length_val = unc_iso['einc_unc']['tof_length']['value']
                self.tof_length_unit = unc_iso['einc_unc']['tof_length']['value_unit']

            self.tof_length_unc = unc_iso['einc_unc']['tof_length']['unc']
            self.tof_length_unc_unit = unc_iso['einc_unc']['tof_length']['unc_unit']
        else:
            self.do_tof_length_unc = False

        if 'tof_t0' in unc_iso['einc_unc']:
            assert hasattr(self, 'tof_length_val'), "Need tof length to propagate tof_err uncertainty"
            self.do_tof_t0_unc = True

            self.tof_t0_val = unc_iso['einc_unc']['tof_t0']['value']
            self.tof_t0_unit = unc_iso['einc_unc']['tof_t0']['value_unit']
            
            if self.tof_t0_val != 0:
                raise ValueError("Setting t0 != 0 currently does not have any impact. Needs implementation")
            
            self.tof_t0_unc = unc_iso['einc_unc']['tof_t0']['unc']
            self.tof_t0_unc_unit = unc_iso['einc_unc']['tof_t0']['unc_unit']
        else:
            self.do_tof_t0_unc = False
        
        # get reference
        self.identifier_iso_deriv1 = unc_iso['einc_unc']['identifier_iso_deriv1']
       #  -------------------------------------------------------------------------------
       
        # loading partial uncertainties in ----------------------------------------------
        self.std_iso = unc_iso['values']
        self.unc_iso_unit = unc_iso['units']
        self.unc_iso_type = unc_iso['type']
        self.unc_iso_typearg = unc_iso['type_arg']
        # -------------------------------------------------------------------------------

        # WIP: error- message if length of Einc =! no_val
      
    # Covariances relative to the cs for main isotope ===================================
    def calculate_cov_iso(self):
        """
        Calculates covariance matrix relative to cs related to unc. of the cs.
        
        output variable:
        cov_iso : covariance matrix relative to cs related to unc. of the cs.
        """
        
        #from CorMatrixShapes import identify_corshape
        import CorMatrixShapes as CMS
        from importlib import reload
        reload(CMS)
        from MatrixFunctions import checking_cov
        from numpy import zeros
        
        cov_iso = zeros([self.no_val,self.no_val],dtype = float)
        sub_cov = zeros([self.no_val,self.no_val],dtype = float)
        sub_cor = zeros([self.no_val,self.no_val],dtype = float)
        
        
        # calculating covariance type --------------------------------------------
        no_unc = 0
        for cor_type in self.unc_iso_type:
            
            # getting cov according to type ......................................
            sub_cor = CMS.identify_corshape(self.no_val,cor_type,self.unc_iso_typearg,no_unc)
            if self.no_val == 1:
                sub_cov[0,0]= self.std_iso[no_unc]*self.std_iso[no_unc]
            else:
                for index1 in range(0,self.no_val):
                    for index2 in range(0,self.no_val):
                        if self.std_iso.ndim > 1:
                            sub_cov[index1,index2] = self.std_iso[index1,no_unc]*self.std_iso[index2,no_unc]*sub_cor[index1,index2]
                        else:
                            sub_cov[index1,index2] = self.std_iso[index1]*self.std_iso[index2]*sub_cor[index1,index2]                        
           # ..................................................................
            
            # testing ..........................................................
            msg = "Testing covariance matrix of type "+str(cor_type)+":"
            print(msg)
            checking_cov(sub_cov)
            print("")
            # ...................................................................
            no_unc+=1        
            cov_iso = cov_iso+sub_cov
        # ------------------------------------------------------------------------
            
        return cov_iso
    # ===========================================================================


# fission cross section absolute data ########################################################################    
class nfcs_absolute(nfcs):
    def __init__(self,general_info,data,unc_iso,features):
        from numpy import sqrt, diag, array
        from ToFCovariance import calculate_cov_eunc_NORATIO
        import WriteOutput as WO
        
        # storing data in appropriate containers -----------------------------------------
        # loading data values, isotope uncertainties and general information in
        super(nfcs_absolute,self).__init__(general_info=general_info,data=data,unc_iso=unc_iso,features=features)
        
        #loading normalization uncertainties in 
        if 'normalization_unc' in unc_iso:
            self.norm_unc_val  = unc_iso['normalization_unc']['value']
            self.norm_unc_unit = unc_iso['normalization_unc']['unit']
        else:
            msg = 'For absolute measurement normalization uncertainty must be provided. If unknown use data as shape data.'
            raise Exception(msg)  
        # --------------------------------------------------------------------------------

        # WIP test if all uncertainties are in %
     
        # calculating the different uncertainty types ------------------------------------        
        # WIP see if all units are the same, i.e. percent
        self.cov_norm = self.calculate_cov_norm()
        self.cov_iso = super(nfcs_absolute,self).calculate_cov_iso()

        # self.cov_eunc = self.calculate_cov_eunc_NORATIO_legacy()
        self.cov_eunc = calculate_cov_eunc_NORATIO(self)

        self.data_total = self.get_data_cor_rel()
        # --------------------------------------------------------------------------------

        # writes output to file ----------------------------------------------------------
        self.data_total['einc_unit'] = self.einc_unit
        self.data_total['values_unit'] = self.values_unit
        self.data_total['einc'] = self.einc           
        print(WO.write_xml_output(self.general_info,self.data_total))
        print(WO.write_json_output_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1))
        
        if self.no_val == 1:
            partial_unc = array([self.einc,self.data_total['data'],self.data_total['rel_unc'][0]*100,sqrt(self.cov_iso[0,0]),\
                                 sqrt(self.cov_norm),sqrt(self.cov_eunc[0,0])])
        else:
            partial_unc = array([self.einc,self.data_total['data'],self.data_total['rel_unc']*100,sqrt(diag(self.cov_iso)),\
                             sqrt(diag(self.cov_norm)),sqrt(diag(self.cov_eunc))])
        partial_unc_header = list(['Einc ('+self.einc_unit+')','Data('+self.values_unit+')', 'Total (%)','Isotope (%)','Normalization unc. (%)','Eunc. Unc. (%)'])
        partial_data = {'header': partial_unc_header,'values':partial_unc}

        print(WO.write_partialunc_output(self.general_info['output_folder'],partial_data))
        print(WO.write_cortofile(self.general_info['output_folder']+'TotalCor.dat',self.einc,self.data_total['cor']))
        # --------------------------------------------------------------------------------
        self.plot_data()

    # covariance due to normalization ====================================================
    def calculate_cov_norm(self):
        """
        Calculates normalization covariance related to unc. of the cs.
        
        output variable:
        cov_norm : covariance matrix relative to cs related to unc. of the cs.
        """
        import CorMatrixShapes as CMS
        from importlib import reload
        reload(CMS)
        from MatrixFunctions import checking_cov
        
        cov_norm = CMS.identify_corshape(self.no_val,'Positive_fully')*self.norm_unc_val*self.norm_unc_val

        
        # testing ..........................................................
        msg = "Testing covariance matrix of type "+str('Positive_fully')+":"
        print(msg)
        checking_cov(cov_norm)
        print("")
        # ...................................................................

        return cov_norm
    # =====================================================================================

  
    # ===========================================================================
    def get_data_cor_rel(self):
        """
        Calculates the total correlation matrix, relative uncertainties and data.
        """
        from MatrixFunctions import checking_cov
        from numpy import zeros, shape, sqrt
        
        data_cor_rel = {}
        
        # get data
        data_cor_rel['data'] = self.values
        
        # get cor ---------------------------------------------------------------
        cov = self.cov_eunc + self.cov_norm + self.cov_iso 
        
         # testing ................................................................
        msg = "Testing total covariance matrix:"
        print(msg)
        checking_cov(cov)
        print("")
        # .......................................................................
                
        dim = shape(cov)[0]
        cor = zeros([dim,dim],dtype = float)
        relunc = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc[index1] = sqrt(cov[index1,index1])/100.0
            for index2 in range(0,dim):
                cor[index1,index2] = cov[index1,index2]/sqrt(cov[index1,index1]*cov[index2,index2])
        
        data_cor_rel['cor'] = cor
        data_cor_rel['rel_unc'] = relunc
        # -----------------------------------------------------------------------              
        return data_cor_rel
    # ===========================================================================

    # Plot shape data ===========================================================
    def plot_data(self):
        """
        Plots the data, uncertainties and covariances.
        """        
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,\
            savefig,errorbar, xscale ,rcParams
        from numpy import diag, sqrt, meshgrid
        from pylab import pcolormesh, colorbar, show, title, imsave
        from MatrixFunctions import get_cor
        
        rcParams.update({'figure.autolayout': True})


        figcount = 0
        # nf cs -----------------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data_log.eps') 
        # -----------------------------------------------------------------------

        # relative uncertainties ------------------------------------------------
        # log x
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k')
            semilogx(self.einc,sqrt(diag(self.cov_iso)),'r')
            semilogx(self.einc,sqrt(diag(self.cov_norm)),'g')
            semilogx(self.einc,sqrt(diag(self.cov_eunc)),'k--')
            legend(['Total','Isotope','Normalization','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        else:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k.')
            semilogx(self.einc,sqrt(self.cov_iso),'r+')
            semilogx(self.einc,sqrt(self.cov_norm),'g>')
            semilogx(self.einc,sqrt(self.cov_eunc),'k<')
        legend(['Total','Isotope','Normalization','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc_log.eps',bbox_inches='tight')
            
        # normal
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            plot(self.einc,self.data_total['rel_unc']*100,'k')
            plot(self.einc,sqrt(diag(self.cov_iso)),'r')
            plot(self.einc,sqrt(diag(self.cov_norm)),'g')
            plot(self.einc,sqrt(diag(self.cov_eunc)),'k')
        else:
            plot(self.einc,self.data_total['rel_unc']*100,'k.')
            plot(self.einc,sqrt(self.cov_iso),'r+')
            plot(self.einc,sqrt(self.cov_norm),'g>')
            plot(self.einc,sqrt(self.cov_eunc),'k<')
        legend(['Total','Isotope','Normalization','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc.eps',bbox_inches='tight')
        # -----------------------------------------------------------------------

        # correlation matrices --------------------------------------------------
        if self.no_val > 1:
            figure(figcount)
            figcount +=1
            x,y = meshgrid(self.einc,self.einc)
            pcolormesh(x,y,self.data_total['cor'],vmin=-1, vmax=1)
            colorbar()
            title('Total Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Total_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_iso)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Isotope Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Isotope_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_norm)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Normalization Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Normalization_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_eunc)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Eunc Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Eunc_Cor.eps')
            show() 
        else:
            print("No correlation matrices plotted as only one data value.")
        
        return None 
    # ===========================================================================
# ###############################################################################
    

# fission cross section absolute data ########################################################################    
class nfcs_shape(nfcs):
    def __init__(self,general_info,data,unc_iso,features):
        from numpy import sqrt, diag, array
        from ToFCovariance import calculate_cov_eunc_NORATIO
        import WriteOutput as WO
        
        # storing data in appropriate containers -----------------------------------------
        # loading data values, isotope uncertainties and general information in
        super(nfcs_shape,self).__init__(general_info=general_info,data=data,unc_iso=unc_iso,features=features)
        
        # checking that no normalization uncertainty is provided -------------------------
        if 'normalization_unc' in unc_iso:
            msg = 'For shape measurements, no normalization and hence associated uncertainties should be provided. If the normalization is known, use data as absolute data.'
            raise Exception(msg)  
        # --------------------------------------------------------------------------------

        # WIP test if all uncertainties are in %
     
        # calculating the different uncertainty types ------------------------------------        
        # WIP see if all units are the same, i.e. percent
        self.cov_iso = super(nfcs_shape,self).calculate_cov_iso()

        # self.cov_eunc = self.calculate_cov_eunc_NORATIO_legacy()
        self.cov_eunc = calculate_cov_eunc_NORATIO(self)

        self.data_total = self.get_data_cor_rel()
        # --------------------------------------------------------------------------------

        # writes output to file ----------------------------------------------------------
        self.data_total['einc_unit'] = self.einc_unit
        self.data_total['values_unit'] = self.values_unit
        self.data_total['einc'] = self.einc           
        print(WO.write_xml_output(self.general_info,self.data_total))
        print(WO.write_json_output_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1))
        
        if self.no_val == 1:
            partial_unc = array([self.einc,self.data_total['data'],self.data_total['rel_unc']*100,sqrt(self.cov_iso),\
                                 sqrt(self.cov_eunc)])
        else:
            partial_unc = array([self.einc,self.data_total['data'],self.data_total['rel_unc']*100,sqrt(diag(self.cov_iso)),\
                             sqrt(diag(self.cov_eunc))])
        partial_unc_header = list(['Einc ('+self.einc_unit+')','Data('+self.values_unit+')', 'Total (%)','Isotope (%)','Eunc. Unc. (%)'])
        partial_data = {'header': partial_unc_header,'values':partial_unc}

        print(WO.write_partialunc_output(self.general_info['output_folder'],partial_data))
        print(WO.write_cortofile(self.general_info['output_folder']+'TotalCor.dat',self.einc,self.data_total['cor']))
        # --------------------------------------------------------------------------------
        self.plot_data()

  
    # ===========================================================================
    def get_data_cor_rel(self):
        """
        Calculates the total correlation matrix, relative uncertainties and data.
        """
        from MatrixFunctions import checking_cov
        from numpy import zeros, shape, sqrt
        
        data_cor_rel = {}
        
        # get data
        data_cor_rel['data'] = self.values
        
        # get cor ---------------------------------------------------------------
        cov = self.cov_eunc +  self.cov_iso 
        
         # testing ................................................................
        msg = "Testing total covariance matrix:"
        print(msg)
        checking_cov(cov)
        print("")
        # .......................................................................
                
        dim = shape(cov)[0]
        cor = zeros([dim,dim],dtype = float)
        relunc = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc[index1] = sqrt(cov[index1,index1])/100.0
            for index2 in range(0,dim):
                cor[index1,index2] = cov[index1,index2]/sqrt(cov[index1,index1]*cov[index2,index2])
        
        data_cor_rel['cor'] = cor
        data_cor_rel['rel_unc'] = relunc
        # -----------------------------------------------------------------------              
        return data_cor_rel
    # ===========================================================================

    # Plot shape data ===========================================================
    def plot_data(self):
        """
        Plots the data, uncertainties and covariances.
        """        
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,\
            savefig,errorbar, xscale ,rcParams
        from numpy import diag, sqrt, meshgrid
        from pylab import pcolormesh, colorbar, show, title, imsave
        from MatrixFunctions import get_cor
        
        rcParams.update({'figure.autolayout': True})


        figcount = 0
        # nf cs -----------------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data_log.eps') 
        # -----------------------------------------------------------------------

        # relative uncertainties ------------------------------------------------
        # log x
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k')
            semilogx(self.einc,sqrt(diag(self.cov_iso)),'r')
            semilogx(self.einc,sqrt(diag(self.cov_eunc)),'g')
            legend(['Total','Isotope','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        else:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k.')
            semilogx(self.einc,sqrt(self.cov_iso),'r+')
            semilogx(self.einc,sqrt(self.cov_eunc),'g>')
            legend(['Total','Isotope','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc_log.eps',bbox_inches='tight')
            
        # normal
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            plot(self.einc,self.data_total['rel_unc']*100,'k')
            plot(self.einc,sqrt(diag(self.cov_iso)),'r')
            plot(self.einc,sqrt(diag(self.cov_eunc)),'g')
        else:
            plot(self.einc,self.data_total['rel_unc']*100,'k.')
            plot(self.einc,sqrt(self.cov_iso),'r+')
            plot(self.einc,sqrt(self.cov_eunc),'g>')
        legend(['Total','Isotope','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc.eps',bbox_inches='tight')
        # -----------------------------------------------------------------------

        # correlation matrices --------------------------------------------------
        if self.no_val > 1:
            figure(figcount)
            figcount +=1
            x,y = meshgrid(self.einc,self.einc)
            pcolormesh(x,y,self.data_total['cor'],vmin=-1, vmax=1)
            colorbar()
            title('Total Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Total_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_iso)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Isotope Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Isotope_Cor.eps')
            show()
            
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_eunc)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Eunc Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Eunc_Cor.eps')
            show() 
        else:
            print("No correlation matrices plotted as only one data value.")
        
        return None 
    # ===========================================================================
# ###############################################################################

# fission cross section absolute data ########################################################################    
class nfcs_cleanratioabsolute(nfcs):
    def __init__(self,general_info,data,unc_iso,reference,features):
        from numpy import sqrt, diag, array, shape, zeros
        from MatrixFunctions import get_cor
        from ToFCovariance import calculate_cov_eunc_RATIO
        import WriteOutput as WO
        
        # storing data in appropriate containers -----------------------------------------
        # loading data values, isotope uncertainties and general information in
        super(nfcs_cleanratioabsolute,self).__init__(general_info=general_info,data=data,unc_iso=unc_iso,features=features)
        
        #loading normalization uncertainties in 
        if 'normalization_unc' in unc_iso:
            self.norm_unc_val  = unc_iso['normalization_unc']['value']
            self.norm_unc_unit = unc_iso['normalization_unc']['unit']
        else:
            msg = 'For absolute clean ratio measurements, a normalization uncertainty must be provided. If unknown use data as shape clean ratio data.'
            raise Exception(msg)  

        # loading reference information in 
        self.reference = reference
        # --------------------------------------------------------------------------------

        # WIP test if all uncertainties are in %
     
        # calculating the different uncertainty types ------------------------------------        
        # WIP see if all units are the same, i.e. percent
        self.cov_norm = self.calculate_cov_norm()
        self.cov_iso = super(nfcs_cleanratioabsolute,self).calculate_cov_iso()

        # self.cov_eunc = calculate_cov_eunc_RATIO_legacy(self)
        self.cov_eunc = calculate_cov_eunc_RATIO(self)

        self.cov_ratio = self.calculate_cov_ratio()
        self.cov_refnum = self.get_cov_refnum()

        self.data_total = self.get_data_cor_rel()
        # --------------------------------------------------------------------------------

        # writes output to file ----------------------------------------------------------
        self.data_total['einc_unit'] = self.einc_unit
        self.data_total['values_unit'] = self.refnum_data_unit
        self.data_total['ratio_unit'] = self.values_unit
        self.data_total['einc'] = self.einc
        self.data_total['ratio'] = self.values

        dim = shape(self.cov_ratio)[0]
        cor_ratio = zeros([dim,dim],dtype = float)
        relunc_ratio = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc_ratio[index1] = sqrt(self.cov_ratio[index1,index1])/100.0
            for index2 in range(0,dim):
                cor_ratio[index1,index2] = self.cov_ratio[index1,index2]/sqrt(self.cov_ratio[index1,index1]*self.cov_ratio[index2,index2])
        
        self.data_total['cor_ratio'] = cor_ratio
        self.data_total['relunc_ratio'] = relunc_ratio
           
        print(WO.write_xml_output(self.general_info,self.data_total))
        print(WO.write_json_output_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1))
        print(WO.write_json_output_ratio_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1,reference))
        
        # For ratio data, not only the total covariance matrix is written to file, but 
        # also the ratio data as the standards folks want to deal with ratio data.
        if self.no_val == 1:
           partial_unc = array([self.einc,self.data_total['data'],self.values,self.data_total['rel_unc']*100,sqrt(self.cov_ratio),sqrt(self.cov_iso),\
                                sqrt(self.cov_refnum),sqrt(self.cov_norm),sqrt(self.cov_eunc)])
        else:
           partial_unc = array([self.einc,self.data_total['data'],self.values,self.data_total['rel_unc']*100,sqrt(diag(self.cov_ratio)),sqrt(diag(self.cov_iso)),\
                             sqrt(diag(self.cov_refnum)),sqrt(diag(self.cov_norm)),sqrt(diag(self.cov_eunc))])
        partial_unc_header = list(['Einc ('+self.einc_unit+')','Data('+self.refnum_data_unit+')', 'Ratio Data('+self.values_unit+')', 'Total (%)','Ratio (%)','Isotope (%)','Reference unc. (%)','Normalization unc. (%)','Eunc. Unc. (%)'])
        partial_data = {'header': partial_unc_header,'values':partial_unc}

        print(WO.write_partialunc_output(self.general_info['output_folder'],partial_data))
        print(WO.write_cortofile(self.general_info['output_folder']+'TotalCor.dat',self.einc,self.data_total['cor']))
        print(WO.write_cortofile(self.general_info['output_folder']+'RatioCor.dat',self.einc,get_cor(self.cov_ratio)))
        # --------------------------------------------------------------------------------
        self.plot_data()

    # covariance due to normalization ====================================================
    def calculate_cov_norm(self):
        """
        Calculates normalization covariance related to unc. of the cs.
        
        output variable:
        cov_norm : covariance matrix relative to cs related to unc. of the cs.
        """
        import CorMatrixShapes as CMS
        from importlib import reload
        reload(CMS)
        from MatrixFunctions import checking_cov
        
        cov_norm = CMS.identify_corshape(self.no_val,'Positive_fully')*self.norm_unc_val*self.norm_unc_val

        
        # testing ..........................................................
        msg = "Testing covariance matrix of type "+str('Positive_fully')+":"
        print(msg)
        checking_cov(cov_norm)
        print("")
        # ...................................................................

        return cov_norm
    # =====================================================================================


      # ===========================================================================
    def get_cov_refnum(self):
        """
        Reads and interpolates covariance matrix of the references nuclear data 
        relative to the reference nuclear data
        
        output variable:
        cov_refnum : covariance matrix of the references nuclear data relative to 
        the reference nuclear data
        """  
        from Manage_ReferenceData import get_reference_datacor
        from MatrixFunctions import checking_cov, get_cor
        from scipy import interpolate
        from numpy import zeros, shape, array, interp, meshgrid, diag, sqrt
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,\
            savefig, xscale, yscale
        from warnings import warn
        from pylab import pcolormesh, colorbar, show, title, imsave
        import BasicPhysicsFunctions as BPF
        from importlib import reload
        reload(BPF)

        # read data and covariance of reference reaction ------------------------
        dataunccor = get_reference_datacor(self.reference)
              
        # testing ................................................................
        if (sum(sum(dataunccor['cor'])) == 0) :
            msg = 'Error: only 0 matrix provided as correlation matrix'
            raise Exception(msg)

        msg = "Testing covariance matrix of reference reaction:"
        print(msg)
        checking_cov(dataunccor['cor'])
        print("")
        # .......................................................................
        # -----------------------------------------------------------------------

        
        
        # convert to adequate lattice -------------------------------------------
        if dataunccor['lattice_unit'] != self.einc_unit:
            dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],self.einc_unit,dataunccor['lattice_unit'])
            dataunccor['lattice_unit'] = self.einc_unit

        if 'latticeunc' in dataunccor:
            lattice_refnum = dataunccor['latticeunc']
            if dataunccor['latticeunc_unit'] != self.einc_unit:
                lattice_refnum = BPF.conversion_to_newUnits(dataunccor['latticeunc'],self.einc_unit,dataunccor['latticeunc_unit'])
                dataunccor['latticeunc_unit'] = self.einc_unit
        else:
            lattice_refnum = dataunccor['lattice']
         
        
        
       # check whether lattice of standard encompasses the exp. lattice ........
        if (min(lattice_refnum) > min(self.einc)):
            msg = "Warning: Minimum outgoing energy of experiment smaller than that of the reference data."
            warn(msg)
        elif (max(lattice_refnum) < max(self.einc)):
            msg = "Warning: Maximum outgoing energy of experiment larger than that of the reference data."
            warn(msg)
        # .......................................................................

        cor = dataunccor['cor']
        relunc = dataunccor['rel_unc']
        dim = shape(cor)[0]
        cov = zeros([dim,dim],dtype = float)
        cov_refnum = zeros([self.no_val,self.no_val],dtype=float)
        self.refnum_data_unit = dataunccor['data_unit']

        for index1 in range(0,dim):
            for index2 in range(0,dim):
                cov[index1,index2] = cor[index1,index2]*relunc[index1]*relunc[index2]
        
        self.refnum_data = array(interp(self.einc,dataunccor['lattice'],dataunccor['data']))


        intcov = interpolate.RectBivariateSpline(lattice_refnum,lattice_refnum,cov, kx=2,ky=2)
        cov_refnum = intcov(self.einc,self.einc)
        # .....................................................................
            
        cor_refnum = get_cor(cov_refnum)
        # -----------------------------------------------------------------------

        # plotting of data for visual tests -------------------------------------

        # n,f ..................................................................
        # as is
        figure(1)
        xscale('log')
        yscale('log')
        semilogx(dataunccor['lattice'],dataunccor['data'],'ks')
        semilogx(self.einc,self.refnum_data,'ro')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel(self.reference['isotope'] + ' (n,f) cs ('+dataunccor['data_unit']+')')
        legend(['Original','Interpolated'])
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation.eps')
        # .......................................................................

        # Rel. unc. .............................................................
        figure(3)
        xscale('log')
        yscale('log')
        semilogx(lattice_refnum,relunc,'ks')
        semilogx(self.einc,sqrt(diag(cov_refnum)),'ro')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel('Relative Uncertainties (%)')
        legend(['Original','Interpolated'])
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation_RelUnc.eps')
        # .......................................................................

        # covariances ...........................................................
        figure(4)
        x,y = meshgrid(self.einc,self.einc)
        pcolormesh(x,y,cor_refnum,vmin=-1, vmax=1)
        colorbar()
        title('Interpolated Cor')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation_Cor.eps')
        show()

        #figure(5)
        #x,y = meshgrid(lattice_refnum,lattice_refnum)
        #pcolormesh(x,y,cor,vmin=-1, vmax=1)
        #colorbar()
        #title('Original Cor')
        #xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        #ylabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        #savefig(self.general_info['output_folder']+'ReferenceDataOrg_Cor.eps')
        #show()
        # .......................................................................

        # -----------------------------------------------------------------------
                
        # check covariance matrix -----------------------------------------------
        msg = "Testing interpolated covariance matrix of reference reaction:"
        print(msg)
        checking_cov(cov_refnum)
        print("")
        # ----------------------------------------------------------------------- 
        
        return cov_refnum
    # ===========================================================================

    
    # ===========================================================================
    def calculate_cov_ratio(self):
        """
        Calculates the covariance matrix of the ratio data.
        """
        from MatrixFunctions import checking_cov
        from numpy import zeros, shape, sqrt
        
        
        # get cov for ratio data --------------------------------------------------
        cov_ratio = self.cov_eunc + self.cov_norm + self.cov_iso 
        
         # testing ................................................................
        msg = "Testing ratio covariance matrix:"
        print(msg)
        checking_cov(cov_ratio)
        print("")
        # .......................................................................
        # -----------------------------------------------------------------------
              
        return cov_ratio
    # ===========================================================================
  
    # ===========================================================================
    def get_data_cor_rel(self):
        """
        Calculates the total correlation matrix, relative uncertainties and data.
        """
        from MatrixFunctions import checking_cov
        from numpy import zeros, shape, sqrt
        
        data_cor_rel = {}
        
        # get data
        data_cor_rel['data'] = self.values*self.refnum_data
        
        # get cor ---------------------------------------------------------------
        cov = self.cov_refnum + self.cov_ratio 
        
         # testing ................................................................
        msg = "Testing total covariance matrix:"
        print(msg)
        checking_cov(cov)
        print("")
        # .......................................................................
                
        dim = shape(cov)[0]
        cor = zeros([dim,dim],dtype = float)
        relunc = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc[index1] = sqrt(cov[index1,index1])/100.0
            for index2 in range(0,dim):
                cor[index1,index2] = cov[index1,index2]/sqrt(cov[index1,index1]*cov[index2,index2])
        
        data_cor_rel['cor'] = cor
        data_cor_rel['rel_unc'] = relunc
        # -----------------------------------------------------------------------              
        return data_cor_rel
    # ===========================================================================

    # Plot shape data ===========================================================
    def plot_data(self):
        """
        Plots the data, uncertainties and covariances.
        """        
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,\
            savefig,errorbar, xscale ,rcParams
        from numpy import diag, sqrt, meshgrid
        from pylab import pcolormesh, colorbar, show, title, imsave
        from MatrixFunctions import get_cor
        
        rcParams.update({'figure.autolayout': True})


        figcount = 0
        # nf cs -----------------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data_log.eps') 
        # -----------------------------------------------------------------------

        # ratio nf cs -----------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.values,yerr=sqrt(diag(self.cov_ratio))*self.values/100.0,fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + '/'+self.reference['isotope']+ ' (n,f) cs')
        savefig(self.general_info['output_folder']+'Ratio.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.values,yerr=sqrt(diag(self.cov_ratio))*self.values/100.0,fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + '/'+self.reference['isotope']+ ' (n,f) cs')
        savefig(self.general_info['output_folder']+'Ratio_log.eps') 
        # -----------------------------------------------------------------------

        # relative uncertainties ------------------------------------------------
        # log x
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k')
            semilogx(self.einc,sqrt(diag(self.cov_iso)),'r')
            semilogx(self.einc,sqrt(diag(self.cov_ratio)),'b')
            semilogx(self.einc,sqrt(diag(self.cov_norm)),'g')
            semilogx(self.einc,sqrt(diag(self.cov_refnum)),'cs')
            semilogx(self.einc,sqrt(diag(self.cov_eunc)),'k--')
            legend(['Total','Isotope','Ratio','Normalization','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        else:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k.')
            semilogx(self.einc,sqrt(self.cov_iso),'r+')
            semilogx(self.einc,sqrt((self.cov_ratio)),'bo')
            semilogx(self.einc,sqrt(self.cov_norm),'g>')
            semilogx(self.einc,sqrt(self.cov_refnum),'cs')
            semilogx(self.einc,sqrt(self.cov_eunc),'k<')
        legend(['Total','Isotope','Ratio','Normalization','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc_log.eps',bbox_inches='tight')
            
        # normal
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            plot(self.einc,self.data_total['rel_unc']*100,'k')
            plot(self.einc,sqrt(diag(self.cov_iso)),'r')
            plot(self.einc,sqrt(diag(self.cov_ratio)),'b')
            plot(self.einc,sqrt(diag(self.cov_norm)),'g')
            plot(self.einc,sqrt(diag(self.cov_refnum)),'cs')
            plot(self.einc,sqrt(diag(self.cov_eunc)),'k')
        else:
            plot(self.einc,self.data_total['rel_unc']*100,'k.')
            plot(self.einc,sqrt(self.cov_iso),'r+')
            plot(self.einc,sqrt((self.cov_ratio)),'bo')
            plot(self.einc,sqrt(self.cov_norm),'g>')
            plot(self.einc,sqrt(self.cov_refnum),'cs')
            plot(self.einc,sqrt(self.cov_eunc),'k<')
        legend(['Total','Isotope','Ratio','Normalization','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc.eps',bbox_inches='tight')
        # -----------------------------------------------------------------------

        # correlation matrices --------------------------------------------------
        if self.no_val > 1:
            figure(figcount)
            figcount +=1
            x,y = meshgrid(self.einc,self.einc)
            pcolormesh(x,y,self.data_total['cor'],vmin=-1, vmax=1)
            colorbar()
            title('Total Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Total_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_iso)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Isotope Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Isotope_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_ratio)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Ratio Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Ratio_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_norm)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Normalization Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Normalization_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_eunc)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Eunc Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Eunc_Cor.eps')
            show() 
        else:
            print("No correlation matrices plotted as only one data value.")
        
        return None 
    # ===========================================================================
#################################################################################

# fission cross section clean ratio shape data ##################################    
class nfcs_cleanratioshape(nfcs):
    def __init__(self,general_info,data,unc_iso,reference,features):
        from numpy import sqrt, diag, array, shape, zeros
        from MatrixFunctions import get_cor
        from ToFCovariance import calculate_cov_eunc_RATIO
        import WriteOutput as WO
        
        # storing data in appropriate containers -----------------------------------------
        # loading data values, isotope uncertainties and general information in
        super(nfcs_cleanratioshape,self).__init__(general_info=general_info,data=data,unc_iso=unc_iso,features=features)
        
        #loading normalization uncertainties in 
        if 'normalization_unc' in unc_iso:
            msg = 'For clean ratio shape measurements, no normalization and hence associated uncertainties should be provided. If the normalization is known, use data as clean ratio absolute data.'
            raise Exception(msg)  

        # loading reference information in 
        self.reference = reference
        # --------------------------------------------------------------------------------

        # WIP test if all uncertainties are in %
     
        # calculating the different uncertainty types ------------------------------------        
        # WIP see if all units are the same, i.e. percent
        self.cov_iso = super(nfcs_cleanratioshape,self).calculate_cov_iso()

        # self.cov_eunc = self.calculate_cov_eunc()
        self.cov_eunc = calculate_cov_eunc_RATIO(self)

        self.cov_ratio = self.calculate_cov_ratio()
        self.cov_refnum = self.get_cov_refnum()

        self.data_total = self.get_data_cor_rel()
        ## --------------------------------------------------------------------------------

        ## writes output to file ----------------------------------------------------------
        self.data_total['einc_unit'] = self.einc_unit
        self.data_total['values_unit'] = self.refnum_data_unit
        self.data_total['ratio_unit'] = self.values_unit
        self.data_total['einc'] = self.einc           
        self.data_total['ratio'] = self.values

        dim = shape(self.cov_ratio)[0]
        cor_ratio = zeros([dim,dim],dtype = float)
        relunc_ratio = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc_ratio[index1] = sqrt(self.cov_ratio[index1,index1])/100.0
            for index2 in range(0,dim):
                cor_ratio[index1,index2] = self.cov_ratio[index1,index2]/sqrt(self.cov_ratio[index1,index1]*self.cov_ratio[index2,index2])
        
        self.data_total['cor_ratio'] = cor_ratio
        self.data_total['relunc_ratio'] = relunc_ratio
 
        print(WO.write_xml_output(self.general_info,self.data_total))
        print(WO.write_json_output_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1))
        print(WO.write_json_output_ratio_EUCLID(self.general_info,self.data_total,self.features,self.identifier_iso_deriv1,reference))
        
        # For ratio data, not only the total covariance matrix is written to file, but 
        # also the ratio data as the standards folks want to deal with ratio data.
        if self.no_val == 1:
           partial_unc = array([self.einc,self.data_total['data'],self.values,self.data_total['rel_unc']*100,sqrt(self.cov_ratio),sqrt(self.cov_iso),\
                                sqrt(self.cov_refnum),sqrt(self.cov_eunc)])
        else:
           partial_unc = array([self.einc,self.data_total['data'],self.values,self.data_total['rel_unc']*100,sqrt(diag(self.cov_ratio)),sqrt(diag(self.cov_iso)),\
                             sqrt(diag(self.cov_refnum)),sqrt(diag(self.cov_eunc))])
        partial_unc_header = list(['Einc ('+self.einc_unit+')','Data('+self.data_total['values_unit']+')', 'Ratio Data('+self.values_unit+')', 'Total (%)','Ratio (%)','Isotope (%)','Reference unc. (%)','Eunc. Unc. (%)'])
        partial_data = {'header': partial_unc_header,'values':partial_unc}

        print(WO.write_partialunc_output(self.general_info['output_folder'],partial_data))
        print(WO.write_cortofile(self.general_info['output_folder']+'TotalCor.dat',self.einc,self.data_total['cor']))
        print(WO.write_cortofile(self.general_info['output_folder']+'RatioCor.dat',self.einc,get_cor(self.cov_ratio)))
        # --------------------------------------------------------------------------------
        self.plot_data()


      # ===========================================================================
    def get_cov_refnum(self):
        """
        Reads and interpolates covariance matrix of the references nuclear data 
        relative to the reference nuclear data
        
        output variable:
        cov_refnum : covariance matrix of the references nuclear data relative to 
        the reference nuclear data
        """  
        from Manage_ReferenceData import get_reference_datacor
        from MatrixFunctions import checking_cov, get_cor
        from scipy import interpolate
        from numpy import zeros, shape, array, interp, meshgrid, diag, sqrt
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,savefig, xscale, yscale
        from warnings import warn
        from pylab import pcolormesh, colorbar, show, title, imsave
        import BasicPhysicsFunctions as BPF
        from importlib import reload
        reload(BPF)

        # read data and covariance of reference reaction ------------------------
        dataunccor = get_reference_datacor(self.reference)
              
        # testing ................................................................
        if (sum(sum(dataunccor['cor'])) == 0) :
            msg = 'Error: only 0 matrix provided as correlation matrix'
            raise Exception(msg)

        msg = "Testing covariance matrix of reference reaction:"
        print(msg)
        checking_cov(dataunccor['cor'])
        print("")
        # .......................................................................
        # -----------------------------------------------------------------------

        
        
        # convert to adequate lattice -------------------------------------------
        if dataunccor['lattice_unit'] != self.einc_unit:
            dataunccor['lattice'] = BPF.conversion_to_newUnits(dataunccor['lattice'],self.einc_unit,dataunccor['lattice_unit'])
            dataunccor['lattice_unit'] = self.einc_unit

        if 'latticeunc' in dataunccor:
            lattice_refnum = dataunccor['latticeunc']
            if dataunccor['latticeunc_unit'] != self.einc_unit:
                lattice_refnum = BPF.conversion_to_newUnits(dataunccor['latticeunc'],self.einc_unit,dataunccor['latticeunc_unit'])
                dataunccor['latticeunc_unit'] = self.einc_unit
        else:
            lattice_refnum = dataunccor['lattice']
         
        
        
       # check whether lattice of standard encompasses the exp. lattice ........
        if (min(lattice_refnum) > min(self.einc)):
            msg = "Warning: Minimum outgoing energy of experiment smaller than that of the reference data."
            warn(msg)
        elif (max(lattice_refnum) < max(self.einc)):
            msg = "Warning: Maximum outgoing energy of experiment larger than that of the reference data."
            warn(msg)
        # .......................................................................

        cor = dataunccor['cor']
        relunc = dataunccor['rel_unc']
        dim = shape(cor)[0]
        cov = zeros([dim,dim],dtype = float)
        cov_refnum = zeros([self.no_val,self.no_val],dtype=float)
        self.refnum_data_unit = dataunccor['data_unit']

        for index1 in range(0,dim):
            for index2 in range(0,dim):
                cov[index1,index2] = cor[index1,index2]*relunc[index1]*relunc[index2]
        
        self.refnum_data = array(interp(self.einc,dataunccor['lattice'],dataunccor['data']))


        intcov = interpolate.RectBivariateSpline(lattice_refnum,lattice_refnum,cov, kx=2,ky=2)
        cov_refnum = intcov(self.einc,self.einc)
        # .....................................................................
            
        cor_refnum = get_cor(cov_refnum)
        # -----------------------------------------------------------------------

        # plotting of data for visual tests -------------------------------------

        # n,f ..................................................................
        # as is
        figure(1)
        xscale('log')
        yscale('log')
        semilogx(dataunccor['lattice'],dataunccor['data'],'ks')
        semilogx(self.einc,self.refnum_data,'ro')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel(self.reference['isotope'] + ' (n,f) cs ('+dataunccor['data_unit']+')')
        legend(['Original','Interpolated'])
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation.eps')
        # .......................................................................

        # Rel. unc. .............................................................
        figure(3)
        xscale('log')
        yscale('log')
        semilogx(lattice_refnum,relunc,'ks')
        semilogx(self.einc,sqrt(diag(cov_refnum)),'ro')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel('Relative Uncertainties (%)')
        legend(['Original','Interpolated'])
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation_RelUnc.eps')
        # .......................................................................

        # covariances ...........................................................
        figure(4)
        x,y = meshgrid(self.einc,self.einc)
        pcolormesh(x,y,cor_refnum,vmin=-1, vmax=1)
        colorbar()
        title('Interpolated Cor')
        xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        ylabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        savefig(self.general_info['output_folder']+'ReferenceDataInterpolation_Cor.eps')
        show()

        #figure(5)
        #x,y = meshgrid(lattice_refnum,lattice_refnum)
        #pcolormesh(x,y,cor,vmin=-1, vmax=1)
        #colorbar()
        #title('Original Cor')
        #xlabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        #ylabel('Incident Neutron Energy ('+dataunccor['lattice_unit']+')')
        #savefig(self.general_info['output_folder']+'ReferenceDataOrg_Cor.eps')
        #show()
        # .......................................................................

        # -----------------------------------------------------------------------
                
        # check covariance matrix -----------------------------------------------
        msg = "Testing interpolated covariance matrix of reference reaction:"
        print(msg)
        checking_cov(cov_refnum)
        print("")
        # ----------------------------------------------------------------------- 
        
        return cov_refnum
    # ===========================================================================

    
    # ===========================================================================
    def calculate_cov_ratio(self):
        """
        Calculates the covariance matrix of the ratio data.
        """
        from MatrixFunctions import checking_cov
        
        
        # get cov for ratio data --------------------------------------------------
        cov_ratio = self.cov_eunc + self.cov_iso 
        
         # testing ................................................................
        msg = "Testing ratio covariance matrix:"
        print(msg)
        checking_cov(cov_ratio)
        print("")
        # .......................................................................
        # -----------------------------------------------------------------------
              
        return cov_ratio
    # ===========================================================================
  
    # ===========================================================================
    def get_data_cor_rel(self):
        """
        Calculates the total correlation matrix, relative uncertainties and data.
        """
        from MatrixFunctions import checking_cov
        from numpy import zeros, shape, sqrt
        
        data_cor_rel = {}
        
        # get data
        data_cor_rel['data'] = self.values*self.refnum_data
        
        # get cor ---------------------------------------------------------------
        cov = self.cov_refnum + self.cov_ratio 
        
         # testing ................................................................
        msg = "Testing total covariance matrix:"
        print(msg)
        checking_cov(cov)
        print("")
        # .......................................................................
                
        dim = shape(cov)[0]
        cor = zeros([dim,dim],dtype = float)
        relunc = zeros(dim,dtype = float)
        for index1 in range(0,dim):
            relunc[index1] = sqrt(cov[index1,index1])/100.0
            for index2 in range(0,dim):
                cor[index1,index2] = cov[index1,index2]/sqrt(cov[index1,index1]*cov[index2,index2])
        
        data_cor_rel['cor'] = cor
        data_cor_rel['rel_unc'] = relunc
        # -----------------------------------------------------------------------              
        return data_cor_rel
    # ===========================================================================

    # Plot shape data ===========================================================
    def plot_data(self):
        """
        Plots the data, uncertainties and covariances.
        """        
        from matplotlib.pyplot import plot,semilogx,legend,xlabel,ylabel,figure,\
            savefig,errorbar, xscale ,rcParams
        from numpy import diag, sqrt, meshgrid
        from pylab import pcolormesh, colorbar, show, title, imsave
        from MatrixFunctions import get_cor
        
        rcParams.update({'figure.autolayout': True})


        figcount = 0
        # nf cs -----------------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.data_total['data'],yerr=self.data_total['data']*self.data_total['rel_unc'],fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + ' Fission Cross Section ('+self.values_unit+')')
        savefig(self.general_info['output_folder']+'Data_log.eps') 
        # -----------------------------------------------------------------------

        # ratio nf cs -----------------------------------------------------------
        # as is
        figure(figcount)
        figcount +=1
        errorbar(self.einc,self.values,yerr=sqrt(diag(self.cov_ratio))*self.values/100.0,fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + '/'+self.reference['isotope']+ ' (n,f) cs')
        savefig(self.general_info['output_folder']+'Ratio.eps') 

        # log x
        figure(figcount)
        figcount +=1
        xscale("log")
        errorbar(self.einc,self.values,yerr=sqrt(diag(self.cov_ratio))*self.values/100.0,fmt='ko')
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel(self.isotope + '/'+self.reference['isotope']+ ' (n,f) cs')
        savefig(self.general_info['output_folder']+'Ratio_log.eps') 
        # -----------------------------------------------------------------------

        # relative uncertainties ------------------------------------------------
        # log x
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k')
            semilogx(self.einc,sqrt(diag(self.cov_iso)),'r')
            semilogx(self.einc,sqrt(diag(self.cov_ratio)),'b')
            semilogx(self.einc,sqrt(diag(self.cov_refnum)),'cs')
            semilogx(self.einc,sqrt(diag(self.cov_eunc)),'k--')
            legend(['Total','Isotope','Ratio','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        else:
            semilogx(self.einc,self.data_total['rel_unc']*100,'k.')
            semilogx(self.einc,sqrt(self.cov_iso),'r+')
            semilogx(self.einc,sqrt((self.cov_ratio)),'bo')
            semilogx(self.einc,sqrt(self.cov_refnum),'cs')
            semilogx(self.einc,sqrt(self.cov_eunc),'k--')
        legend(['Total','Isotope','Ratio','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc_log.eps',bbox_inches='tight')
            
        # normal
        figure(figcount)
        figcount +=1
        if self.no_val > 1:
            plot(self.einc,self.data_total['rel_unc']*100,'k')
            plot(self.einc,sqrt(diag(self.cov_iso)),'r')
            plot(self.einc,sqrt(diag(self.cov_ratio)),'b')
            plot(self.einc,sqrt(diag(self.cov_refnum)),'cs')
            plot(self.einc,sqrt(diag(self.cov_eunc)),'k--')
        else:
            plot(self.einc,self.data_total['rel_unc']*100,'k.')
            plot(self.einc,sqrt(self.cov_iso),'r+')
            plot(self.einc,sqrt((self.cov_ratio)),'bo')
            plot(self.einc,sqrt(self.cov_refnum),'cs')
            plot(self.einc,sqrt(self.cov_eunc),'k--')
        legend(['Total','Isotope','Ratio','Reference','Eunc'],\
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
        xlabel('Incident Neutron Energy ('+self.einc_unit+')')
        ylabel('Relative Uncertainties (%)')
        savefig(self.general_info['output_folder']+'RelUnc.eps',bbox_inches='tight')
        # -----------------------------------------------------------------------

        # correlation matrices --------------------------------------------------
        if self.no_val > 1:
            figure(figcount)
            figcount +=1
            x,y = meshgrid(self.einc,self.einc)
            pcolormesh(x,y,self.data_total['cor'],vmin=-1, vmax=1)
            colorbar()
            title('Total Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Total_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_iso)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Isotope Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Isotope_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_ratio)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Ratio Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Ratio_Cor.eps')
            show()
            
            figure(figcount)
            figcount +=1
            cor = get_cor(self.cov_eunc)
            pcolormesh(x,y,cor,vmin=-1, vmax=1)
            colorbar()
            title('Eunc Cor')
            xlabel('Incident Neutron Energy ('+self.einc_unit+')')
            ylabel('Incident Neutron Energy ('+self.einc_unit+')')
            savefig(self.general_info['output_folder']+'Eunc_Cor.eps')
            show() 
        else:
            print("No correlation matrices plotted as only one data value.")
        
        return None 
    # ===========================================================================



#################################################################################

# query for uncertainty of reference without normalization

# fission cross section absolute data ########################################################################    
#class nfcs_indirectratioabsolute(nfcs):
# query for uncertainties of reference including normalization, maybe, two trsls and energy unc. could be given

# fission cross section absolute data ########################################################################    
#class nfcs_indirectratioshape(nfcs):

# query for uncertainty of reference without normalization
