from numpy import array, concatenate, searchsorted, trapz, sqrt, log, zeros_like, linspace, exp
from scipy.interpolate import interp1d

def get_integrated_model(edges, interpolator):
    """Calculates a vector of integrals of model the defined by bin edges and a function (interpolator).

    Args:
        edges (array-like): Bin edges on x-axis.
        interpolator (scipy.interpolate.interp1d): interpolator object that has attribute interpolator.x and returns y values when called as interpolator(x).

    Returns:
        array-like: Vecotr of model integrated values over bins, output likely needs to be normalized by diff(edges)
    """
    # Compute integral over each bin using trapezoidal rule
    binned_model = []

    for i in range(len(edges) - 1):
        a, b = min(edges[i:i+2]), max(edges[i:i+2])

        # assume interpolator has fine samples within edges for accuracy
        index_bounds = searchsorted(interpolator.x, a), searchsorted(interpolator.x, b)
        index_bounds = [min(index_bounds), max(index_bounds)]

        x_fine = concatenate([array([a]), interpolator.x[index_bounds[0]:index_bounds[1]], array([b])])
        
        y_fine = interpolator(x_fine)
        integral = trapz(y_fine, x_fine)

        binned_model.append(integral)

    return array(binned_model)


def fwhm_to_sigma(fwhm):
    """Convert FWHM to standard deviation for a Gaussian."""
    return fwhm / (2 * sqrt(2 * log(2)))

def broaden_with_energy_resolution_vector(
    E_fine, sigma_fine, 
    E_exp, FWHM_exp, 
    kernel_extent=5.0, # how many sigmas to cover
    n_kernel_points=1001 # integration points per kernel
):
    """
    Parameters:
      E_fine:        1D array of theoretical energies (MeV), fine grid
      sigma_fine:    1D array of theoretical cross sections
      E_exp:         1D array of experimental energies at which broadened model will be caluclated
      FWHM_exp:      1D array of experimental energy resolutions (FWHM, MeV)
      kernel_extent: How many sigmas to cover in kernel integration
      n_kernel_points: Number of integration points in kernel
    Returns:
      sigma_broad:   1D array, same length as E_exp
    """
    # Create an interpolator for the theoretical cross section
    sigma_interp = interp1d(E_fine, sigma_fine, kind='linear', bounds_error=False, fill_value=0.0)
    sigma_broad = zeros_like(E_exp)
    
    for k, (E0, FWHM) in enumerate(zip(E_exp, FWHM_exp)):
        sigma_G = fwhm_to_sigma(FWHM)
        # Integration range: [E0 - N*sigma, E0 + N*sigma]
        E_min = E0 - kernel_extent * sigma_G
        E_max = E0 + kernel_extent * sigma_G
        E_kernel = linspace(E_min, E_max, n_kernel_points)
        # Gaussian kernel, normalized
        kernel = exp(-0.5*((E_kernel - E0)/sigma_G)**2)
        kernel /= trapz(kernel, E_kernel)  # normalize kernel to 1
        
        # Evaluate cross section at these points
        sigma_vals = sigma_interp(E_kernel)
        
        # Numerically integrate
        sigma_broad[k] = trapz(sigma_vals * kernel, E_kernel)
    
    return sigma_broad