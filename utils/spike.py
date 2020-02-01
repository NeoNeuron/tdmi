# author : Kyle Chen
# created : 9-Sep-19
"""
Modules for processing spike trains in neural data

"""
import numpy as np

def spike(raster, index = 0, tmin = 0, tmax = None, dt = 0.5, savefile = None):
    """
    Convert spike train to binary sequence and output to *.npy file.
    
    Parameters
    ----------
    raster : narray
        Original raster data
    index : int
        Index of target neuron
    tmin : float, optional
        Lower bound of time range, default 0
    tmax : float, optional
        Upper bound of time range,
        default: the maximum time point in raster
    dt : float
        Binning size, default 0.5 ms
    savefile : string, optional
        Name of output file.
        If none, then no saving to files

    Returns
    -------
    binary_spike : ndarray of bools
        The array of binary spiking sequence.

    Example
    -------
    >>> raster_dat = np.array([[0,0.4],[1,2.3],[0,5.0],[0,7.6]])
    ...  spike(raster_dat, index = 0, dt = 1)
    array([1, 0, 0, 0, 0, 1, 0, 1])

    """
    # get single neuron spike train
    single_spike = raster[raster[:,0]==index,1]
    if tmax == None:
        tmax = np.floor(single_spike[-1] + dt)

    # index-lize the time range
    t_range = np.array([tmin,tmax])
    t_range = np.floor(t_range/dt).astype(int)
    
    # create and process binary sequence
    binary_spike = np.full(t_range[1]-t_range[0], False, dtype=np.int32)
    binary_id = np.floor(single_spike/dt).astype(int)
    # truncate the binary sequence
    binary_id = binary_id[np.all([binary_id>=t_range[0], binary_id<=t_range[1]], axis=0)]
    binary_spike[binary_id-t_range[0]] = True

    # save file
    if savefile != None:
        np.save(savefile, binary_spike)
    return binary_spike

def main(prefix, raster_fname, spike_fname, id, dt, time_range, verbose):
    '''
    Generate single spike train from raster data file.

    Parameters
    ----------
    prefix : string
        prefix of input and output data
    raster_fname : string
        filename of raster data
    spike_fname : string
        filename of spike train data
    id : int
        id of neuron
    dt : float
        sampling time step
    time_range : list of floats, with len(time_range)=2
        time range of single spike train
    verbose : bool
        flag to print the runing time

    Return
    ------
    single_spike : ndarray of int32
        single spike train

    '''
    import time
    t_start = time.time()
    raster = np.genfromtxt(prefix + raster_fname, delimiter=',')
    single_spike = spike(raster, index = id, tmin = time_range[0], tmax = time_range[1], dt = dt, savefile = prefix+spike_fname)
    t_finish = time.time()
    if verbose:
        print('>> Output spike train:')
        print(single_spike)
        print('>> Calculating TDMI took %5.3e s'% (t_finish-t_start))
    return single_spike

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "Single spike train generator.")
    parser.add_argument('prefix', type = str, help = 'working directory of source data and output data')
    parser.add_argument('raster_fname', type = str, help = 'filename of raster data, *.csv')
    parser.add_argument('spike_fname', type = str, help = 'filename of spike train data, *.npy')
    parser.add_argument('--id', type = int, default=0, help = 'id of neuron (default: 0)')
    parser.add_argument('--dt', type = float, default=0.5, help = 'time step of sampling (default: 0.5)')
    parser.add_argument('-t', '--time_range', type = float, nargs=2, default=[0, None], help = 'time range of spike train (default: [0, None])')
    parser.add_argument('-v', '--verbose', action='store_true', help = 'enable verbose of running time')
    args = parser.parse_args()
    print(args)
    main(**(vars(args)))
