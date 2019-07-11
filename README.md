# Stochastic-Volatility-Models

##### R Code to accompany 

####  _An Approach to Efficient Fitting of Univariate and Multivariate Stochastic Volatility Models_

#### The paper is available here: https://www.stat.pitt.edu/stoffer/xxx.pdf

For the bibTeX item for this site, I used the following:
<pre>
@misc{GitGongStoffer2019,
  author = {Gong, Chen and Stoffer, David S.},
  title = {{Stochastic Volatility Models}},
  howpublished = "\url{https://github.com/nickpoison/Stochastic-Volatility-Models/}",
  year = {2019}, 
  note = "[GitHub Repository]"
}  
</pre>



* The data are in the folder *data* and are compressed R data files.
* The various PGAS files are in the folder *R* ... these are sourced in the files used to run the examples.
* Each example is identified by starting with `run_` and then a self describing title.  You just run the code, it will call the data file and PGAS procedure as needed.

-------------------

 You'll need the following R packages to run all the code:

* astsa
* plyr
* MASS 
* mcmc   