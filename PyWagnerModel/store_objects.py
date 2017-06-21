import pickle
import os, errno

def silent_remove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def store_objects(obj, filename):
    with open(filename, 'ab') as output:
        pickle.dump(obj, output, -1)
