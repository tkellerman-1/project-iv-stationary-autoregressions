The code files listed here are used to generate posterior densities for the autoregressive model of order p, or AR(p) model, and the time-varying equivalent, TVAR(p) model.
This includes two Gibbs samplers for the AR(p) model, one of which contains a Metropolis sampling step. There are also stan files to validate both of these, as well as Stan files to compute the posterior and forecast distributions for the TVAR(p) model.
The Metropolis sampling scheme in the file PACMetropolis2 (1).R contains contributions from Dr Sarah Heaps at the University of Durham.
