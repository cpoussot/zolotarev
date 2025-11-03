# Loewner Framework for the Zolotarev Z3 and Z4 problems

## Overview

This page accompany a the paper entitled "The Loewner Framework applied ro Zolotarev signand ratio problems", by A.C. Antoulas, I.V. Gosea and C. Poussot-Vassal. In this work, we propose a numerical study concerning the (rational) approximation of functions connected to the 3rd and 4th Zolotarev problems. We compare numerical results for various methods, including the Loewner framework (LF), but also the standard AAA algorithm and a recently proposed extensions of AAA (namely sign and Lawson). We show that the Loewner framework is fast and reliable, and provides approximants with a high level of accuracy, sometimes even more accurate than near-optimal ones for higher degrees. Last but not least, the Loewner framework is a direct method, for which the running time is significantly lower than that of the iterative AAA-Lawson method. Moreover, for the latter, the running time increases substantially with the degree of the approximant, whereas for the Loewner method, it does not. We also show that on specific topologies, the LF allows to recover the poles/zeros distriubution with a machine precision.

This page provides the necessary numerical tools allowing to reproduce the results of the paper.

## Contribution claim

The contributions and highlights of this paper, discovered through a comprehensive numerical comparison performed here, are listed below:
- we show that the LF solves Zolotarev problems by compressing the number of interpolation points. This approach does not need any iterations and yields solutions quite close to the optimal ones;
- we conduct an extensive numerical study and comparison of the performance of different methods (w.r.t. to computing time, accuracy of fit and interpretability of the polynomials, for approximating several sign functions defined on various domains showing, among others, the computational advantage of LF;
- we emphasize the LF, and demonstrate that this method yields, in most cases, very accurate results in a fast and reliable way, with no user intervention and no iteration. As a matter of fact, we highlight the LF as a valid alternative to these methods;
- we show that the rational functions computed with LF come very close to optimal approximants and, in some specific cases (notably the symmetric ones), recovers the structure of the optimal solution as well as the symmetry property, i.e., the eigenvalues and zeros distributions, whereas the other methods tend to add spurious and badly distributed poles and zeros. Such a feature is indeed decisive to discover the property of the topology when only data values are available;
- we point out that in simplified cases for which the exact solution is known, the LF recovers the polynomial structure (i.e., realness, alternation of coefficients in polynomial representation), whereas AAA adds complex, non-trivial artifacts that can’t be easily explained.

![mLF](doc/spiral2_3D.png "Z4 for spiral")

![mLF](doc/pm2_2D.png "Z4 for Pac Man")


## Main reference

```
@article{AGPV:2025,
	Author	= {A.C. Antoulas and I.V. Gosea and C. Poussot-Vassal},
	Doi 	= {},
	Journal = {submitted},
	Number 	= {},
	Pages 	= {},
	Title 	= {{The Loewner Framework applied ro Zolotarev signand ratio problems}},
	Volume 	= {},
  	Month   = {},
	Year 	= {},
	Note    = {}, 
}
```

# The "zol" MATLAB package 

The code (`+zol` folder)  provided in this GitHub page is given for open science purpose. Its principal objective is to accompany the readers, and thus aims at being as educative as possible rather than industry-oriented. Evolutions (numerical improvements) may come with time. Please, cite the reference above if used in your work and do not hesitate to contact us in case of bug of problem when using it. Below we present an example of use, then functions list are given.
It is also meant to allow reproduction of the results in the paper.

## Simple MATLAB code examples

We provide a series of simple codes that describe how to deploy the LF and how to compare with some AAA approaches. These demo files are meant to reproduce some results contained in the above mentionned paper. More specifically, we include
- `demo0_LF.m`: runs the LF to solve Z3 and Z4 problems (line 10 select the case name, line 11 provide the desired order). 
- `demo1_LF_vs_AAA_1ab.m`: runs the LF and AAA to  solve Z3 and Z4 problems `1a` (two circles) and `1b` (two real lines). Note that here you need to add the `chebfun` package (line 10); the latter may be downloaded at the [Chebfun main page](https://www.chebfun.org/). As for these two cases the optimal solution is known, attention is pushed to the numerators, denominators, poles and zeros obtained by each methods. As pointed in the present paper, we highlight the acurate results in term of monomial and poles/zeros reached by the LF.
- `demo2_LF_vs_AAA.m`: compares performances and poles/zeros the LF and AAA to  solve Z3 and Z4 the collection of proposed problems. Here attention is given to the accuracy, poles and zeros and computational time.
- `demo2_LF_vs_AAA_time.m`: compares performances and poles/zeros the LF and AAA to solve Z3 and Z4 the collection of proposed problems. Here, attention is given to the ratio number and computational time for each methods.


## Functions description

Please check help in the functions below, where all informations and details are provided: 
```Matlab
help zol.example
help zol.example2data
help zol.loewner
help zol.pb4_to_pb3
```

## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.


## Disclaimer

Once again, this deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced code, to deal with practical and theoretical issues/limitations is currently under developpement by the authors.





