# spectral_clustering
This is an implementation of the normalized spectral clustering algorithm.<br/>
A description of the algorithm can be found in [this paper](https://ai.stanford.edu/~ang/papers/nips01-spectral.pdf).<br/>
<br/>
The code is written in both Python & C, which are integrated using the Python/C API.<br/>
<br/>

## install & run
- Download the repository to your computer.<br/>
- Install the C extension by running ``` python setup.py build_ext --inplace ```<br/>
- Run ``` python spkmeans.py [#-of-clusters] spk [features.csv]```<br/>
