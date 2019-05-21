def convert_numpy_bytes_to_unicode(runtimes):
    '''
    Our old runtime arrays use the S dtype descriptor, which
    in Python 3 becomes bytes, making str comparison
    more cumbersome.  We can just change the dtype of
    the string fields to be Unicode.  Doing so on
    a structured arrays is also cumbersome, hence
    this function.

    Our new runtime arrays just use Unicode by default.
    '''
    runtimes_descr = runtimes.dtype.descr
    runtimes_descr = [list(d) for d in runtimes_descr]
    for d in runtimes_descr:
        if 'S' in d[1]:
            d[1] = d[1].replace('S','U')
    runtimes_descr = [tuple(d) for d in runtimes_descr]
    runtimes = runtimes.astype(runtimes_descr)

    return runtimes
