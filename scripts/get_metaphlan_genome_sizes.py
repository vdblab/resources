import pickle
import bz2
import sys

if __name__ == "__main__":
    try:
        sys.stderr.write("parsing pickle\n")
        with bz2.open(sys.argv[1], 'rb') as f:
            data = pickle.load(f)
        sys.stderr.write("writing items\n")
        for k, v in data["taxonomy"].items():
            sys.stdout.write(f"{k}\t{v[0]}\t{v[1]}\n") 
    except Exception as e:
        raise(e)
