# Stochastic-Volatility-Models

---
## [The univariate case is included in `astsa` as `SV.mcmc`](https://github.com/nickpoison/astsa/blob/master/README.md)
---

<br/><br/>
<br/><br/>

 R Code to accompany the Sept 2020 and final version of

  _A Note on Efficient Fitting of Stochastic Volatility Models_

The paper has been published online: [jtsa.12561](https://onlinelibrary.wiley.com/doi/10.1111/jtsa.12561)


* The data are in the folder *data* and are compressed R data files.
* The various PGAS files are in the folder *R* ... these are sourced in the files used to run the examples.
* Each example is identified by starting with `run_` and then a self describing title.  You just run the code, it will call the data file and PGAS procedure as needed.
* Added an example from stochvol R package, but the essential part
of the code is in one of the vignettes.

 _You'll need the following R packages to run all the code:_

* astsa
* plyr
* MASS 
* mcmc  
* stochvol (_needed only to run their example, figure 2_)

<br/> 



----
----
The bibTeX entry for the current version is:

```
@article{doi:10.1111/jtsa.12561,
author = {Gong, Chen and Stoffer, David S.},
title = {A Note on Efficient Fitting of Stochastic Volatility Models},
journal = {Journal of Time Series Analysis},
year = {2021},
volume = {42},
number = {2},
pages = {186-200},
keywords = {Ancestral sampling, efficient Markov chain Monte Carlo, particle Gibbs, stochastic volatility},
doi = {10.1111/jtsa.12561},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12561},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/jtsa.12561},
}
```

And plain text is

```
Gong, C. and Stoffer, D.S. (2021), A Note on Efficient Fitting of Stochastic Volatility Models. 
          J. Time Ser. Anal., 42: 186-200. https://doi.org/10.1111/jtsa.12561
```


---

For the bibTeX item to the code here, I used the following:

```
@misc{GitGongStoffer2020,
author = {Gong, Chen and Stoffer, David S.},
title = {{Stochastic Volatility Models}},
howpublished = "\url{https://github.com/nickpoison/Stochastic-Volatility-Models/}",
month = {09},
year = {2020}, 
note = "[GitHub Repository]"
}  
```
