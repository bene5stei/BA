def __store_as_pickle(obj, filename):

    import pickle
    from os.path import isdir

    ofile = open(filename, 'wb')
    pickle.dump(obj, ofile)

    if isdir(filename):
        print(f"created: {filename}")
