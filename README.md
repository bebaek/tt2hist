# tt2hist

Calculate a 2-, 3-, or 4-dimensional correlation function from time-tagged event data. Originally used to calculate g(2), g(3), and g(4) of photon statistics: ["High-order temporal coherences of chaotic and laser light"](https://doi.org/10.1364/OE.18.001430) and more.

- Input datafile format is proprietary so such header and cpp files are not included here. However, it should be possible to adapt the code for an open time-tagged detection data format.

- Tested with MS Visual C++ Express 2010. Both 32- and 64-bit work.

- As in many other scientific programs, this was written in C++ for a performance reason. Parallel computing could be the next step to improve this further.
