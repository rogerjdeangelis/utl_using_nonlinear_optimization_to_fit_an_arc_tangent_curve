Using Nonlinear Optimization to fit an arc tangent curve

y = 2*a*tan(x/a)/pi

Two Solutions

    1. WPS/proc Python and SAS proc nlin
    2. WPS/Proc R   (Optimizes the squared error f(x)-observed Y Grid search?)

see
https://tinyurl.com/y8l9r8wq
https://github.com/rogerjdeangelis/utl_using_nonlinear_optimization_to_fit_an_arc_tangent_curve

There are other methods to fit the simple curve but this demonstrates an
WPS/proc R/Python  and SAS/Stat nlin solutions.
I realize the derivative is not needed but wanted to include for more
reasonable examples.

Framing a good question is much harder than the answer.

see
https://stackoverflow.com/questions/50189130/approximating-a-curve-with-arctan

Weihuang Wong profile
https://stackoverflow.com/users/6455166/weihuang-wong



INPUT  (Find the best value of a to fit the equation below)
=====

y = 2*atan(x/a)/pi

WORK.HAVE total obs=50

Obs      X        Y

  1     1000    0.101
  2     2000    0.169
  3     3000    0.209
  4     4000    0.256
  5     5000    0.289
...
 47    47000    0.994
 48    48000    0.991
 49    49000    0.995
 50    50000    0.993


     |
1.00 +                      ******
     |                 ******
     |              ***
     |             ***
     |            *
0.75 +           *
     |          *
     |         * Optimum value for a
  Y  |        **
     |       **  WPS proc python and SAS NLIN   a=7069.4
0.50 +       *   WPS/Proc R                     a=7069.4
     |
     |      *
     |     *
     |     *
0.25 +    *
     |    *
     |   *
     |   *
     |
0.00 +
     ---+---------+---------+---------+--
        0       20000     40000     60000

                       X

Example Output
                          Approx         Approximate 95%
Parameter  Estimate    Std Error        Confidence Limits

a            7069.4        391.7      6282.4      7856.5




PROCESS
=======

WPS/PROC PYTHON AND SAS NLIN
----------------------------

We use Python to differentiate y = 2*a*tan(x/a)/pi with respect to a.
We will use the derivative in SAS Proc NLIN

%utl_submit_wps64("
options set=PYTHONHOME 'C:\Users\backup\AppData\Local\Programs\Python\Python35\';
options set=PYTHONPATH 'C:\Users\backup\AppData\Local\Programs\Python\Python35\lib\';
proc python;
submit;
from sympy import *;
import pandas as pd;
from sympy import *;
x, a = symbols('x,a', real=True);
dif=diff(2*atan(x/a)/pi, a);
print(dif);
print(solve(dif, x));
endsubmit;
run;quit;
");

Derivative

-2*x/(pi*a**2*(1 + x**2/a**2))

proc nlin best=10 plot method=marquardt data=have;
  parms a=500 to 20000 by 50;
  model y=2*atan(x/a)/constant('pi');
  der.a=-2*x/(constant('pi')*a**2*(1 + x**2/a**2));
  output out=b p=yhat r=yresid;
options ls=64 ps=32;

* rename to get plot with width less than 64;
proc plot data=b(rename=
(y=y1234567890123456789012345
yresid=y12345678901234567890)
);
  plot y1234567890123456789012345*x='*' yhat*x='p' /overlay vpos=25;
  plot y12345678901234567890*x / vref=0 vpos=25;
run;quit;

options ls=171 ps=66;

/* output

                              Approx         Approximate 95%
Parameter Estimate    Std Error        Confidence Limits

a           7069.4        391.7      6282.4      7856.5
*/


WPS PROC R   WPS/PROC PYTHON AND SAS NLIN
===========

%utl_submit_wps64('
libname sd1 "d:/sd1";
options set=R_HOME "C:/Program Files/R/R-3.3.1";
libname wrk sas7bdat "%sysfunc(pathname(work))";
proc r;
submit;
source("C:/Program Files/R/R-3.3.1/etc/Rprofile.site", echo=T);
library(haven);
have<-read_sas("d:/sd1/have.sas7bdat");
head(have);
xs<-have$X;
ys<-have$Y;
form <- function(a, xs) 2 * atan(xs / a) / pi;
loss <- function(a, xs, data) mean((form(a, xs) - ys)^2);
optimized <- optim(1000, loss, xs = xs, data = data, method = "Brent", lower = 1, upper = 10000);
a <- optimized$par;
a;
endsubmit;
run;quit;
');

/*
The WPS System

[1] 7069.42
*/

OUTPUT
======

WPS/PROC PYTHON AND SAS NLIN
----------------------------

       Estimation Summary

Method                  Marquardt
Iterations                      4
R                        2.082E-6
PPC(a)                   8.073E-7
RPC(a)                   6.342E-6
Object                   2.33E-10
Objective                0.278885
Observations Read              50
Observations Used              50
Observations Missing            0

NOTE: An intercept was not specified for this model.

                              Sum of      Mean           Approx
Source                  DF   Squares    Square  F Value  Pr > F

Model                    1   32.6748   32.6748  5740.95  <.0001
Error                   49    0.2789   0.00569
Uncorrected Total       50   32.9536


                              Approx         Approximate 95%
Parameter      Estimate    Std Error        Confidence Limits

a                7069.4        391.7      6282.4      7856.5



                              Approx
Parameter      Estimate    Std Error    Approximate 95% Confidence Limits

a                7069.4        391.7      6282.4      7856.5




1.00 +                      ******
     |    FIT          *****
     |                ***   pppp
     |              **  ppppp
     |             *pppp
     |            *pp
0.75 +          p*
     |         p*
     |        p*
2345 |       pp
     |       p**
     |      p*
0.50 +      p
     |     p *
     |      *
     |     **
     |    p
     |     *
0.25 +    *
     |    *
     |   *
     |
     |   *
     |
0.00 +
     ---+---------+---------+---------+--
        0       20000     40000     60000

                       X


890 |
    |
    |
    |    RESIDUALS
    |
0.1 +                   A   B A
    |                 BAABBB BABB
    |              AAB A
    |             BAA
    |            A
    |   A
0.0 +---A-------BA----------------------
    |          A
    |          A
    |    A   AA
    |    A
    |     A AAA
0.1 +     AAA
    |      A
    |
    |
    |
    |
0.2 +
    |
    ---+---------+---------+---------+--
       0       20000     40000     60000

                 X


WPS/PROC R
----------

The WPS System

[1] 7069.42

*                _               _       _
 _ __ ___   __ _| | _____     __| | __ _| |_ __ _
| '_ ` _ \ / _` | |/ / _ \   / _` |/ _` | __/ _` |
| | | | | | (_| |   <  __/  | (_| | (_| | || (_| |
|_| |_| |_|\__,_|_|\_\___|   \__,_|\__,_|\__\__,_|

;
xs <- seq(1000,50000,by=1000);
data <- c(0.101, 0.169, 0.209, 0.256, 0.289, 0.373, 0.391, 0.418,
0.477, 0.528, 0.579, 0.570, 0.602, 0.657, 0.690, 0.720, 0.747,
0.764, 0.781, 0.811, 0.842, 0.847, 0.864, 0.889, 0.871, 0.894,
0.897, 0.915, 0.926, 0.931, 0.932, 0.940, 0.950, 0.963, 0.956,
0.967, 0.967, 0.963, 0.971, 0.980, 0.986, 0.988, 0.986, 0.985,
0.993, 0.992, 0.994, 0.991, 0.995, 0.993);

form <- function(a, xs) 2 * atan(xs / a) / pi
loss <- function(a, xs, data) mean((form(a, xs) - data)^2)
optimized <- optim(1000, loss, xs = xs, data = data, method = "Brent", lower = 1, upper = 10000)
a <- optimized$par

*          _       _   _
 ___  ___ | |_   _| |_(_) ___  _ __  ___
/ __|/ _ \| | | | | __| |/ _ \| '_ \/ __|
\__ \ (_) | | |_| | |_| | (_) | | | \__ \
|___/\___/|_|\__,_|\__|_|\___/|_| |_|___/

;

see process

