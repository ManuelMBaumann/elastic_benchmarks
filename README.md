[![DOI](https://zenodo.org/badge/64120949.svg)](https://zenodo.org/badge/latestdoi/64120949)

# Two benchmark problems for the time-harmonic elastic wave equation in 2D and 3D

We present two benchmark problems for the time-harmonic elastic wave equation in 2D and 3D. The aim of this repository is to make the test cases considered in Section 5 of [Baumann et al., 2016](http://www.ewi.tudelft.nl/en/the-faculty/departments/applied-mathematics/reports/) publicly available.   

## Mathematical formulation

We consider the elastic wave equation in a frequency-domain formulation,

![elastic wave eqn](/figs/main_eqn.jpg)

where the unknown **u** is the displacement vector at the k-th frequency. The code we present allows an inhomogeneous medium, and implements a stress-free material-air boundary condition in the north, and first-order Sommerfeld radiation boundary conditions elsewhere.

A finite-element discretization is described in all detail in Section 2 of [Baumann et al., 2016](http://www.ewi.tudelft.nl/en/the-faculty/departments/applied-mathematics/reports/). The discrete problem yields,

![Alt text](/figs/discr_eqn.jpg)

In contrast to the work in our publication, the resulting linear systems are here solved sequentially (over k) using python's sparse solver `scipy.sparse.linalg.spsolve` without preconditioning. For large problems, we recommend to set the flag `storing=True` which stores the discretization matrices **K,C,M** in a folder `/matrix_io` and does *not* solve the resulting linear system. The matrices are then stored in [matrix market](http://math.nist.gov/MatrixMarket/) format, with a MATLAB interface `matrix_io/matlab_io.m` provided.

![Alt text](/figs/spy_plot.png)

The above spy plots can be obtained by setting `spy=True`. The left plot resembles a problem of `ndims=2` and finite element `degree=2`. The right plot belongs to a problems of `ndims=3`and `degree=1`. 

## A sample of numerical test cases

### 2D elastic wedge problem
We define a new *elastic* wedge problem consisitng of three layers with varying physical parameters given in a computational domain **[0,600]x[0,1000]** meters. The parameters are given in the following table, and Lame parameters are computed [accordingly](http://scienceworld.wolfram.com/physics/LameConstants.html).

|parameter    | layer #1 | layer #2 | layer #3 | 
|-------------|----------|----------|----------|
|rho [kg/m^3] | 1800     | 2100     | 1950     |
|c_p [m/s]    | 2000     | 3000     | 2300     |
|c_s [m/s]    | 800      | 1600     | 1100     |

The python script `elast_wedge.py` with the following parameters,

`python3 elast_wedge.py --ndims=2 --freq=[16.0] --dx=2.5 --dz=2.5 --plots=True --nprocs=4`

yields the numerical results:

![Alt text](/figs/wedge_2d_plots.png)

### 3D elastic wedge problem
We extend the 2D wedge problem to 3D by expending all parameters in y-direction. A 3D discretization can be obtained by setting the flag `ndims=3` and by specifying `dy`.

`python3 elast_wedge.py --ndims=3 --freq=[4.0] --dx=10.0 --dy=10.0 --dz=10.0 --storing=True --nprocs=4`

![Alt text](/figs/wedge_3d_plots.png)

### Marmousi-II problem
In the folder `/data` the data set of the [Marmousi-II](https://www.google.nl/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&ved=0ahUKEwiw2OqZlofOAhUUM8AKHdzNAcoQFggvMAM&url=http%3A%2F%2Fmcee.ou.edu%2Faaspi%2Fpublications%2F2006%2Fmartin_etal_TLE2006.pdf&usg=AFQjCNFIhgpermjt0pQo7m51uuHtnWrqIg&sig2=hQNDPOGvwUyhq85xi0nH_g&cad=rja) problem is stored in [segy](https://docs.obspy.org/master/packages/obspy.io.segy.html) format. We consider the following subset of the problem:

![Alt text](/figs/marm2_cs.png)

Solving the Marmousi-II at `freq=6` Hertz can be done via:

`python3 marmousi2.py --freq=[6.0] --n_course=4 --nprocs=4`

Here, we use every fourth grid point of the original problem in each spatial direction which yields `dx=dz=5.0`. The source is located at **(Lx/3,0)**.

![Alt text](/figs/marm2_f6.png)

## Usage and installation

This code is purely written in [python3](https://www.python.org/download/releases/3.0/) using standard numerical libraries such as NumPy and SciPy. For the finite element discretization we use [nutils](http://www.nutils.org/).

The following installation steps are necessary:

* Clone this repository: `git clone https://github.com/ManuelMBaumann/elastic_benchmarks.git`
* Install [nutils](http://www.nutils.org/) via `pip install git+https://github.com/joostvanzwieten/nutils@955bc67d219496e26b037f47855709a222850d7c`
* Run your first, low-frequency, 2D test case `python3 elast_wedge.py --ndims=2 --dx=10.0 --dz=10.0 --freq=[4.0] --nprocs=4`, and view your results with `firefox ~/public_html/elast_wedge.py/latest/log.html`.

For the Marmousi-II problem, two additional steps are required:

* Download the Marmousi-II data set using the provided script `download_marmousi2.sh` [~ 450 Mb].
* Install [obspy](https://github.com/obspy/obspy/wiki) via `pip install obspy`. In particular, we use the segy reader `obspy.io.segy.core._read_segy` for loading the data sets. Obspy requires [libxml2](http://xmlsoft.org/downloads.html) and [libxslt](https://git.gnome.org/browse/libxslt/).


A more detailed description on the installation of nutils can be found in this [document](http://joostvanzwieten.github.io/nutils-by-example/).

**Note:** All plots will be saved as .png in a folder `~/public_html/elast_wedge.py/latest/`. In the same folder, a file `log.html` contains information about the program evaluation and embeds all figures.

## Declaration

The [author](http://www.manuelbaumann.de) is a PhD student in Numerical Analysis at TU Delft. My research is focused on linear solvers and preconditioning. I highly encourage experts in geophysics to comment on the numerical results and to [get in touch](mailto:m.m.baumann@tudelft.nl).

## References

* [M. Baumann, R. Astudillo, Y. Qiu, E. Ang, M.B. Van Gijzen, and R-E. Plessix. A Preconditioned Matrix Equation Approach for the Time-Harmonic Elastic Wave Equation at Multiple Frequencies. TU Delft Technical Report (2016).](/literature/msss_idr_elast_report.pdf)
* [M. Baumann and M.B. Van Gijzen. Nested Krylov methods for shifted linear systems. SIAM J. Sci. Comput., 37(5), S90-S112 (2015).](/literature/140979927-2.pdf)



