# PhaseReconstruction
This is a class project I did for UCSC's AM214 Applied Dynamical Systems. In this repository are a few Jupyter notebooks that walk through numerical examples using Taken's theory of phase space reconstruction on various dynamical systems. 

## What is Phase Space Reconstruction?
Phase space reconstruction is a method of analyzing chaotic dynamical systems that involves simplifying the data at hand. For systems where we have data across several variables, but little or no knowledge of the equations governing the system, it can be challenging to extract meaningful information from the large amount of data. Phase space reconstruction involves taking a single dimension of the system and creating a series of vectors using a time delayed sampling. These vectors are then used to analyze the system. 

The first step for performing phase space reconstruction is creating the time delayed vectors. This involves choosing a dimension of the system to sample from, a dimension, and a time delay. The dimension is how many elements are in each vector, and the time delay (often represented as tau) is how much time passes between each element in a given vector. It's worth noting that some dimensions and time delays are better than others, so in practice its best to try a few different values. Say that we've chosen a dimension to sample that has one data point each second, and that we choose dimension of 3 and tau of 5. Then, our vectors would be constructed as (t=0, t=5, t=10), (t=1, t=6, t=11), (t=2, t=7, t=12), ... , (t, t+tau, t+2tau). 

After creating the vectors, we can plot them and analyze any structure that is present. As an example, I've pulled a figure from the `example_lorenz.ipynb` notebook. The first figure below is the phase space diagram of a lorenz oscillator, with parameters such that the two oscillatory lobes are present. 

Observe that the second image -- a plot after reconstruction has been done -- still maintains the same structure as the original phase space. 

## How to use this code
This repository is intended to provide examples that illustrate phase space reconstruction. If you do run them, also play with different values of tau and the dimensionality and see what different structures you can uncover. \

The reconstruction itself is all done in the `reconstructPhase()` function. In addition to specifying a dimension and time delay tau, the function has a `density` parameter which can be used to throw out some of the reconstructed vectors. Doing so doesn't change the end result, it just makes things run a litle faster. 

Feel free to modify or distribute this code as you see fit. 
