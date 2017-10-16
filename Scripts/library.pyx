# distutils: language = c++

import cython
cimport numpy as np
import numpy

ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE_64_t
ctypedef np.int32_t DTYPE_int_t
ctypedef np.int64_t DTYPE_int64_t
ctypedef np.uint32_t DTYPE_uint_t
ctypedef np.int8_t DTYPE_int8_t
cdef double Inf = numpy.inf
cdef double NaN = numpy.nan

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double log2(double x) nogil
    double log10(double x) nogil
    double sqrt(double x) nogil
    double pow(double x, double x) nogil
    double abs(double x) nogil
    double round(double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binning_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins not None,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_t, ndim=1] signal not None,
        double chrom_mean,
        int startfend,
        int startfend1,
        int stopfend1,
        int startfend2,
        int stopfend2,
        int startbin1,
        int stopbin1,
        int startbin2,
        int stopbin2,
        int diag):
    cdef long long int fend1, fend2, afend1, afend2, j, k, index, map1, map2
    cdef long long int start1, start2, stop1, stop2
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    cdef int diag2 = diag * 2
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5)) - diag
    cdef int num_parameters = fend_indices.shape[1]
    with nogil:
        for fend1 in range(startfend1, stopfend1):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            k = 0
            index = map1 * (num_bins - 1) - map1 * (map1 + 1 - diag2) / 2 - 1 + diag
            if map1 == startbin1:
                start2 = max(startfend2, fend1 + 2)
            else:
                start2 = fend1 + 2
            if map1 == stopbin1:
                stop2 = stopfend2
            else:
                stop2 = num_fends
            for fend2 in range(start2, stop2):
                map2 = mapping[fend2]
                if map2 < 0 or (diag == 0 and map2 == map1):
                    continue
                if fend2 - fend1 == 3 and fend1 % 2 == 0:
                    continue
                 # give starting expected value
                value = 1.0
                # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                for j in range(num_parameters):
                    afend1 = fend1 + startfend
                    afend2 = fend2 + startfend
                    if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                        value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                    else:
                        value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                # if finding enrichment, correct for distance
                if not parameters is None:
                    distance = log(<double>(mids[fend2] - mids[fend1]))
                    while distance > parameters[k, 0]:
                        k += 1
                    value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                signal[index + map2] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binning_expected2(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins not None,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices not None,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] signal not None,
        np.ndarray[DTYPE_int_t, ndim=2] ranges not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices1 not None,
        int startfend,
        double chrom_mean):
    cdef long long int fend1, fend2, afend1, afend2, i, j, k, bin1, bin2, map1, map2
    cdef double distance, value
    cdef int num_parameters = fend_indices.shape[1]
    cdef long long int num_bins = indices0.shape[0]
    with nogil:
        for i in range(num_bins):
            bin1 = indices0[i]
            bin2 = indices1[i]
            for fend1 in range(ranges[bin1, 0], ranges[bin1, 1]):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                k = 0
                for fend2 in range(max(fend1 + 2, ranges[bin2, 0]), ranges[bin2, 1]):
                    map2 = mapping[fend2]
                    if map2 == -1:
                        continue
                    if fend2 - fend1 == 3 and (fend1 + startfend) % 2 == 0:
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    for j in range(num_parameters):
                        afend1 = fend1 + startfend
                        afend2 = fend2 + startfend
                        if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                            value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                        else:
                            value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding enrichment, correct for distance
                    if not parameters is None:
                        distance = log(<double>(mids[fend2] - mids[fend1]))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    signal[i] += value
    return None



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_complexity_variance(
        np.ndarray[DTYPE_64_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=2] valid not None,
        np.ndarray[DTYPE_64_t, ndim=2] var not None,
        int width):
    cdef long long int i, j, x, y, N
    cdef double mean, mean2, temp
    cdef long long int num_bins = data.shape[0]
    cdef long long int window = var.shape[1]
    with nogil:
        for i in range(num_bins - width * 2 + 1):
            for j in range(i + width, min(i + window, num_bins) - width + 1):
                N = 0
                mean = 0.0
                mean2 = 0.0
                for x in range(i, i + width):
                    for y in range(j, j + width):
                        if valid[x, y] == 1:
                            N += 1
                            temp = data[x, y]
                            mean += temp
                            mean2 += temp * temp
                if N > 1:
                    mean /= N
                    mean2 /= N
                    var[i, j - i] = mean2 - mean * mean
                else:
                    var[i, j - i] = NaN
    return None
