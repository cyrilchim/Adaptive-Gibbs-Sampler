## Description

This manual describes how to write a C++ target density functions using template.hpp in order to run `AMCMC` function.

Armadillo (RcppArmadillo) library was extensively used to write the package: http://arma.sourceforge.net


### Global variables

The following global variables are defined in [Adaptive_Gibbs.hpp](../include/Adaptive_Gibbs.hpp) and are available for users needs:

  - `int dim` - dimensionality of the target distribution;

  - `int par` - number of parameters/blocks to update by the Gibbs/MwG algorithms. If no blocking is used, par = d;

  - `int burn_in` - number of burn-in steps;

  - `double lowerb` - lower bound <a href="https://www.codecogs.com/eqnedit.php?latex=$&space;\epsilon&space;$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$&space;\epsilon&space;$" title="$ \epsilon $" /></a> from the the ARSG algorithm. By defualt, it is set to be `0.1/par^2`;

  - `double am` - sequence <a href="https://www.codecogs.com/eqnedit.php?latex=$&space;a_m&space;$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$&space;a_m&space;$" title="$ a_m $" /></a> from the ARSG algorithm;

 - `double m_iter` - index `m` of the sequence <a href="https://www.codecogs.com/eqnedit.php?latex=$&space;a_m&space;$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$&space;a_m&space;$" title="$ a_m $" /></a>;

  - `uvec blocking_structure` - `par+1` dimensional vector that allows easy access to the elements of blocks for the blocking parameter of AMCMC(...) function. Defined as `blocking_structure[0] = 0, blocking_structure[i] - blocking_structure[i-1] = blocking[i-1], i=0 ,.., par`;

  - `vec current_location` - the current location of the chain;

  - `string working_directory` - write/read directory, where the input of the chain is stored. You can set the directory using `set_working_direcoty(directory)` function described in [AMCMC_appendix.md](../man/AMCMC_appendix.md)

  - `int err_type` - error variable to be used in the constructor of inherited classes of distribution_class. See gaussian_target.hpp for an example of use.  Non-zero value exits the program and throws the corepsonding error `err_name`;

  - `string err_name` - error reference name;

## C++ target density

The target density is defined as an inherited class of distribution_class. The main functions which may or may not be specified depending on the parameters of [AMCMC(...)](../man/AMCMC_info.md) are:

  - `double density(vec theta)`  - target density computed at `theta`;

  - `double logdensity(vec theta)` - logarithm of the target density at `theta`;
  
  - `double full_cond(vec theta, int block_ind)` - full conditional of the block with index `block_ind`;
  
  - `double logfull_cond(vec theta, int block_ind)`  - logarithm of the full conditional of the block wih index `block_ind`;
 
  - `vec sample_full_cond(vec theta, int block_ind)` - sample block `block_ind` from its full conditional distribution

Note, that the densities should be specified up to a normalising constant. At least one of the functions needs to be specified in order to run the algorithm


### Example: gaussian target

As an example we specify all the functions of the `distribution_class` for the multivariate gaussian distribution. It is strongly recommended to read [gaussian_target.hpp](../examples/gaussian_target.hpp) file in order to be able to write your own target density. 

Notice the flexibility of the classes. The user is free to specify constructor `gaussian()` and `~destructor()` to list a set of operations done before and after running the sampling algorithms. For example, a set of additional parameters is introduced in the body of the inherited class: 

```C++
mat S, I, Q, A, D; // S - covariance matrix of the target. Q - precision matrix. I - identity matrix
field<mat> A_block; // A set of matrices such that `A_block[i,j] = (I - diag(Q_{11},..,Q_{ss}) Q)_{ij}`, where `s` is the number of blocks in the Gibbs sampling scheme;
field<mat> sd, inv_sd; // set of matrices such that sd(i) = Q_{ii}^{-1} and inv_sd = sd^{-1};
```

We followed here the notations from `Roberts and Sahu (1997)`. All these parameters were specified in the constructor `gaussian::gaussian()`.

The variable `err_type` was used in the constructor `gaussian::gaussian()`. It precludes a situation of using a non-existent or wrongly defined covariance matrix (which has to be specified through `set_covariance(matrix)` or `set_example_covariance(dim)` functions defined in gaussian_target.hpp)

Note the ease of defining the target density:
```C++
double gaussian::density(vec theta)
{
    return exp(-1./2*dot(theta.t()*Q,theta));
    
}
```


## C++ user defined target density

1. To access the elements of a block gp , one can refer to the corresponding subvector of theta through blocking_structure parameter:
```C++
theta.subvec(blocking_structure(block_ind), blocking_structure(block_ind+1)-1)
```
2. In order to specify sampling from the full conditionals function `sample_full_cond(vec theta, int block_ind)`, R/Rcpp random number generator should be used. For example, to sample a 100 standard normal random variables type `rnorm(100)`, which is of NumericVector type. Use `runif(100)` to sample the uniform random variables.

3. The user is expected to be able to mimic gaussian_target.hpp file and specify the necessary functions in template.hpp.

To call  AMCMC(...) function for the user defined C++ density, run
```C++
AMCMC(distribution_type = "new_density", ... )
```

In order to change a name of the target distribution density to, say, `my_distribution`:

1. Create a new .hpp file and call it, say, `my_distribution.hpp`.
2. Copy paste the content of template.hpp to the newly created file `my_distribution.hpp`. 
3. Change the name of the class from `new_density` to `my_distribution`.
4. Open file density_list.hpp In the preamble add 
 ```C++
 #include "my_distribution .hpp"
 ```
5. Open file density_list.hpp and add

```C++
 else if(distribution_type == "my_distribution")
    {
      c_density = new my_distribution;
    }
```

Now you can call 
```C++
AMCMC(distribution_type = "my_distribution", ... )
```
