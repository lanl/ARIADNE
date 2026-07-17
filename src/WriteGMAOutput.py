# write a function that writes GMA output


def WriteGMAOutput(GMAnewnumber,year,authors,journal,general_info,uncertainty_descriptor,data,reference,\
                   abs_shape,abc_components,Eunc,Enrsl,uncertainties,controlnumber,normalizationunc):
    
    from numpy import array, shape, ones, zeros, loadtxt, interp,arange,sqrt, delete,max, size
    import numbers

    # open GMA-outputfile --------------------------------------------------------------------
    filename = general_info['output_folder']+'DS'+GMAnewnumber+'.CRD'
    f=open(filename,'w')
    # ----------------------------------------------------------------------------------------
    
    # Write first line -----------------------------------------------------------------------
    [reaction_name,reaction_no,reaction_ids] = translate_reaction_to_GMA_code(general_info,reference,abs_shape)
    #f.write('{0:>4s}{0:4s}{0:24s}'.format((GMAnewnumber,year,reaction_name)))
    f.write('{0:>4}'.format(GMAnewnumber.rjust(4))); f.write('{0:4s}'.format(year)); f.write('{0:24s}'.format(reaction_name.ljust(24)))
    f.write('{0:28}'.format(authors.ljust(28)));f.write('{0:10}'.format(journal.ljust(10))); f.write('\n')

    # ----------------------------------------------------------------------------------------
    
    # Write second line ----------------------------------------------------------------------
    f.write('{0:2s}'.format('1'.rjust(2)));f.write('{0:2d}'.format(reaction_no))
    f.write('{0:2s}'.format('0'.rjust(2))); # i.e., always give partial uncertainties
    # 1 would mean give total uncertainties and correlation matrix.
    f.write('{0:2s}'.format('0'.rjust(2))); # i.e., no correlation to other data sets
    # EXTEND LATER!!!!!!!!!!!!
    f.write('{0:3d}'.format(shape(uncertainty_descriptor)[0]))
    
    if size(data['einc']) > 1:
        f.write('{0:5d}'.format(shape(data['einc'])[0]))
    else:
        f.write('{0:5d}'.format(1))
    for reaction_id in reaction_ids:
        f.write('{0:3d}'.format(reaction_id))
    f.write('\n')
    # ----------------------------------------------------------------------------------------
    
    # Write comments out ---------------------------------------------------------------------
    for comment in uncertainty_descriptor:
        f.write(' ');f.write(comment);f.write('\n')    
    # ----------------------------------------------------------------------------------------
    
    # write normalization unc. -> only for absolute data -------------------------------------
    if abs_shape == 'absolute':
        for norm in normalizationunc:
            print(f.write('{0:5.1f}'.format(norm)))
        
        tags = zeros(10,dtype=int)
        for index in arange(shape(normalizationunc)[0]):
            if normalizationunc[index] != 0.0:
                tags[index] = 1
        for tag in tags:
            f.write('{0:3d}'.format(tag))
        f.write('\n')
    # ----------------------------------------------------------------------------------------

    # write a,b,c components -----------------------------------------------------------------
    for line in abc_components:
        for abc in line:
            f.write('{0:5.2f}'.format(abc))
        f.write('\n')   
    # ----------------------------------------------------------------------------------------
    
    # write identifier line -> not used ------------------------------------------------------
    # control number 9 should always be used for total uncertainties or statistical uncertainties
    # control number 0 means the uncertainty source is zero.
    if controlnumber[2] != 9:
        print("WARNING: Should the third control number be 9?")
        
    for index in arange(0,11):
        f.write('{0:3d}'.format(controlnumber[index]))
    f.write('\n')
    # ----------------------------------------------------------------------------------------
    
    # write data and partial uncertainties ---------------------------------------------------
    if size(data['einc']) > 1:
        for index in arange(shape(data['einc'])[0]):
            #f.write('{0:10.4g}'.format(data['einc'][index]))
            f.write(f'{data['einc'][index]:10.4g}')
            #f.write('{0:10.4f}'.format(data['values'][index]))
            f.write(f'{data['values'][index]:10.4g}')
            f.write('{0:5.1f}'.format(Eunc[index]));f.write('{0:5.1f}'.format(Enrsl[index]))
            for unc in uncertainties[index,:]:
                f.write('{0:5.1f}'.format(unc))
            f.write('\n')
    else:
        for index in arange(1):
            f.write('{0:10.4f}'.format(data['einc']));f.write('{0:10.4f}'.format(data['values']))
            f.write('{0:5.1f}'.format(Eunc))
            f.write('{0:5.1f}'.format(Enrsl[index]))
            for unc in uncertainties[index,:]:
                f.write('{0:5.1f}'.format(unc))
            f.write('\n')
    # ----------------------------------------------------------------------------------------
    
    # end line -------------------------------------------------------------------------------
    f.write('EBEBEBEBEBEB'); f.write('\n')
    # ----------------------------------------------------------------------------------------
    
    f.close()
    return filename

def translate_reaction_to_GMA_code(general_info,reference,abs_shape):
    # This function translate ARIADNE isotope specification to GMA specifications.
    reaction_ids = [0,0,0,0,0]
    
    # defining reaction number as of Table 2 of ANL-NDM-139. ----------------------------------
    if reference == 'none':
        if abs_shape == 'absolute':
            reaction_no = 1
        elif abs_shape == 'shape':
            reaction_no = 2    
        reaction_ids[0] = matching_isotope_name(general_info['isotope'])[1]
    else:
        if abs_shape == 'absolute':
            reaction_no = 3
        elif abs_shape == 'shape':
            reaction_no = 4
        reaction_ids[0] = matching_isotope_name(general_info['isotope'])[1]
        reaction_ids[1] = matching_isotope_name(reference['isotope'])[1]
    # -----------------------------------------------------------------------------------------
    
    # defining reaction_name ------------------------------------------------------------------
    if reference == 'none':
        reaction_name = matching_isotope_name(general_info['isotope'])[0]+'('+general_info['reaction']+')'
    else:
        reaction_name = matching_isotope_name(general_info['isotope'])[0]+'('+general_info['reaction']+')/'+\
        matching_isotope_name(reference['isotope'])[0]+'('+reference['reaction']+')'
        
    # -----------------------------------------------------------------------------------------
       
    return [reaction_name,reaction_no,reaction_ids]

def matching_isotope_name(isotopein):
    isotopeout = 'not_found'

    if isotopein == 'Pu-239':
        isotopeout = 'Pu9'
        isotopeid = 9
    elif isotopein == 'U-235':
        isotopeout = 'U5'
        isotopeid = 8
    elif isotopein == 'U-238':
        isotopeout = 'U8'
        isotopeid = 10
    elif isotopein == "Li-6":
        isotopeout = 'LI6'
        isotopeid = 1 
    elif isotopein == "B-10":
        isotopeout = 'B10'
        isotopeid = 3
    return isotopeout,isotopeid
