
This project is built upon the work in [SALT](https://github.com/chengengjie/salt), which was created and open-sourced by [chengengjie](https://github.com/chengengjie).
I would like to thank the original authors for their open-source contributions. You can find the original code here: [SALT](https://github.com/chengengjie/salt)。

# Delay-Driven Rectilinear Steiner Tree Construction

This project aims to develop a rectilinear Steiner tree construction algorithm designed to reduce Elmore delay meanwhile maintaining a bounded wirelength.

## Dependencies

* g++ (version >= 5.4.0) or other working c++ compliers
* CMake (version >= 3.5.1)
* Boost (version >= 1.58)
* Python (version 3, optional, for utility scripts)

## Quick Start

The simplest way to build and run is as follows.
~~~
$ mkdir build
$ cd build
$ cmake ..
$ make -j16
$ cd ../run
~~~

The executable file is located in the `run` directory. 

You can draw the results of the tree using `plot/draw.py`. An example result:

![salt](/toys/DelayTree_toy1.tree.png)

## Benck
The nets extracted from [ICCAD'15 Contest Problem B](https://doi.org/10.1109/ICCAD.2015.7372672) can be downloaded via [Dropbox](https://www.dropbox.com/sh/gcq1dh84ko9rjpz/AAAVT0pLZG_FMiOi0ORiKddva?dl=0). 
You can randomly assign slacks to pins by running the executable file `run/benchmark_constuct`.


## Note

* Make sure you have files `POWV9.dat` and `POST9.dat` in the executable's directory, if not copy them from path `src/salt/base/flute`.
* The source of a net should be guaranteed to have an index of 0. 