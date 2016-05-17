# runDM

*With runDMC, It's Tricky. With runDM, it's not.*

`runDM` is a tool for calculating the running of the couplings of Dark Matter (DM) to the Standard Model (SM) in simplified models with vector mediators. By specifying the mass of the mediator and the couplings of the mediator to SM fields at high energy, the code can be used to calculate the couplings at low energy, taking into account the mixing of all dimension-6 operators. The code can also be used to extract the operator coefficients relevant for direct detection, namely low energy couplings to up, down and strange quarks and to protons and neutrons. Further details about the physics behind the code can be found in Appendix B of [arXiv:1605.XXXXX](http://arxiv.org/abs/1605.XXXXX).

At present, the code is written in two languages: *Mathematica* and *Python*. If you are interested in an implementation in another language, please get in touch and we'll do what we can to add it. But if you want it in Fortran, you better be ready to offer something in return. Installation instructions and documentation for the code can be found in `doc/manual.pdf`. We also provide a number of example files:

- For the Python code, we provide an example script as well as Jupyter Notebook. A static version of the notebook can be viewed [here](http://nbviewer.jupyter.org/github/bradkav/runDM/blob/master/python/runDM-examples.ipynb).

- For the Mathematica code, we provide an example notebook. We also provide an example of how to interface with the [NRopsDD code](http://www.marcocirelli.net/NROpsDD.html) for obtaining limits on general models.

If you make use of `runDM` in your own work, please cite it as:

>F. D’Eramo, B. J. Kavanagh & P. Panci (2016). runDM (Version X.X) [Computer software]. Available at https://github.com/bradkav/runDM/

making sure to include the correct version number. Please also cite the associated papers:

>A. Crivellin, F. D’Eramo & M. Procura, New Constraints on Dark Matter Effective Theories from Standard Model Loops, Phys. Rev. Lett. 112 (2014) 191304 [arXiv:1402.1173],

>F. D’Eramo & M. Procura, Connecting Dark Matter UV Complete Models to Direct Detection Rates via Effective Field Theory, JHEP 1504 (2015) 054 [arXiv:1411.3342],

>F. D’Eramo, B. J. Kavanagh & P. Panci, You can hide but you have to run: direct detection with vector mediators, (2016) [arXiv:1605.XXXXX].

Please contact Bradley Kavanagh (bradkav@gmail.com) for any questions, problems, bugs and suggestions.
