Springbok

I wrote this module to model collective cell migration. This code is posted for
the sake of other researchers. Because I was the only one running this code, it
is not as thoroughly documented as it would be if I were actively collaborating
with others. Therefore, feel free to get in touch to get it up and running.

The code consists of four modules:
* springbok models collective cell motion
* tiger solves the diffusion equation
* platypus is a Matplotlib interface
* meerkat coordinates parameter studies
All four modules are under the MIT license and available on Github.

springbok/springbok models a general ensemble of cells.
`springbok/n_signal_relay` has a variety of scripts that set the parameters for
the different studies in the MBoC paper. Critically,
`springbok/n_signal_relay/n_exo.py` has the code that determines how the
neutrophils sense and affect the distribution of signaling molecules.

Please consider citing the article, "Determining whether observed eukaryotic
cell migration indicates chemotactic responsiveness or random chemokinetic
motion," available here,
https://sites.google.com/site/alexszatmary/articles
which establishes how to relate gradient conditions and cell directionality.
