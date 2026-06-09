import argparse
import numpy as np
from io import StringIO

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--direction', type=int, default=1, help="flipping direction (0: x, 1: y, 2: z)")
parser.add_argument('-i', '--input', type=str, required=True, help='bvec file to flip')
parser.add_argument('-o', '--output', type=str, required=True, help='bvec file name to save')
args = parser.parse_args()
matrix = np.loadtxt(args.input, dtype='f') #, delimiter='  ')
matrix[args.direction,:]=-matrix[args.direction,:]
np.savetxt(args.output,matrix,fmt='%.8f')