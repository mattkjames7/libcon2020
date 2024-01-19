import numpy as np
from scipy.special import j0,j1

def main():
    """
    Save scipy bessel function output to compare to C++ output.
    """
    x = np.linspace(0,100,1001).astype("float64")
    z0 = j0(x).astype("float64")
    z1 = j1(x).astype("float64")

    with open("scipy-bessel.bin","w") as f:
        np.int32(x.size).tofile(f)
        x.tofile(f)
        z0.tofile(f)
        z1.tofile(f)


if __name__ == "__main__":
    main()