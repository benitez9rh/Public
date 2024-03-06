# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:35:47 2022

@author: s2132627
"""



""" Distribution check """


'''
Distributions option are:
    
alpha: An alpha continuous random variable.
anglit: An anglit continuous random variable.
arcsine: An arcsine continuous random variable.
argus: Argus distribution
beta: A beta continuous random variable.
betaprime: A beta prime continuous random variable.
bradford: A Bradford continuous random variable.
burr: A Burr (Type III) continuous random variable.
burr12: A Burr (Type XII) continuous random variable.
cauchy: A Cauchy continuous random variable.
chi: A chi continuous random variable.
chi2: A chi-squared continuous random variable.
cosine: A cosine continuous random variable.
crystalball: Crystalball distribution
dgamma: A double gamma continuous random variable.
dweibull: A double Weibull continuous random variable.
erlang: An Erlang continuous random variable.
expon: An exponential continuous random variable.
exponnorm: An exponentially modified Normal continuous random variable.
exponweib: An exponentiated Weibull continuous random variable.
exponpow: An exponential power continuous random variable.
f: An F continuous random variable.
fatiguelife: A fatigue-life (Birnbaum-Saunders) continuous random variable.
fisk: A Fisk continuous random variable.
foldcauchy: A folded Cauchy continuous random variable.
foldnorm: A folded normal continuous random variable.
genlogistic: A generalized logistic continuous random variable.
gennorm: A generalized normal continuous random variable.
genpareto: A generalized Pareto continuous random variable.
genexpon: A generalized exponential continuous random variable.
genextreme: A generalized extreme value continuous random variable.
gausshyper: A Gauss hypergeometric continuous random variable.
gamma: A gamma continuous random variable.
gengamma: A generalized gamma continuous random variable.
genhalflogistic: A generalized half-logistic continuous random variable.
genhyperbolic: A generalized hyperbolic continuous random variable.
geninvgauss: A Generalized Inverse Gaussian continuous random variable.
gilbrat: A Gilbrat continuous random variable.
gompertz: A Gompertz (or truncated Gumbel) continuous random variable.
gumbel_r: A right-skewed Gumbel continuous random variable.
gumbel_l: A left-skewed Gumbel continuous random variable.
halfcauchy: A Half-Cauchy continuous random variable.
halflogistic: A half-logistic continuous random variable.
halfnorm: A half-normal continuous random variable.
halfgennorm: The upper half of a generalized normal continuous random variable.
hypsecant: A hyperbolic secant continuous random variable.
invgamma: An inverted gamma continuous random variable.
invgauss: An inverse Gaussian continuous random variable.
invweibull: An inverted Weibull continuous random variable.
johnsonsb: A Johnson SB continuous random variable.
johnsonsu: A Johnson SU continuous random variable.
kappa4: Kappa 4 parameter distribution.
kappa3: Kappa 3 parameter distribution.
ksone: Kolmogorov-Smirnov one-sided test statistic distribution.
kstwo: Kolmogorov-Smirnov two-sided test statistic distribution.
kstwobign: Limiting distribution of scaled Kolmogorov-Smirnov two-sided test statistic.
laplace: A Laplace continuous random variable.
laplace_asymmetric: An asymmetric Laplace continuous random variable.
levy: A Levy continuous random variable.
levy_l: A left-skewed Levy continuous random variable.
levy_stable: A Levy-stable continuous random variable.
logistic: A logistic (or Sech-squared) continuous random variable.
loggamma: A log gamma continuous random variable.
loglaplace: A log-Laplace continuous random variable.
lognorm: A lognormal continuous random variable.
loguniform: A loguniform or reciprocal continuous random variable.
lomax: A Lomax (Pareto of the second kind) continuous random variable.
maxwell: A Maxwell continuous random variable.
mielke: A Mielke Beta-Kappa / Dagum continuous random variable.
moyal: A Moyal continuous random variable.
nakagami: A Nakagami continuous random variable.
ncx2: A non-central chi-squared continuous random variable.
ncf: A non-central F distribution continuous random variable.
nct: A non-central Student’s t continuous random variable.
norm: A normal continuous random variable.
norminvgauss: A Normal Inverse Gaussian continuous random variable.
pareto: A Pareto continuous random variable.
pearson3: A pearson type III continuous random variable.
powerlaw: A power-function continuous random variable.
powerlognorm: A power log-normal continuous random variable.
powernorm: A power normal continuous random variable.
rdist: An R-distributed (symmetric beta) continuous random variable.
rayleigh: A Rayleigh continuous random variable.
rice: A Rice continuous random variable.
recipinvgauss: A reciprocal inverse Gaussian continuous random variable.
semicircular: A semicircular continuous random variable.
skewcauchy: A skewed Cauchy random variable.
skewnorm: A skew-normal random variable.
studentized_range: A studentized range continuous random variable.
t: A Student’s t continuous random variable.
trapezoid: A trapezoidal continuous random variable.
triang: A triangular continuous random variable.
truncexpon: A truncated exponential continuous random variable.
truncnorm: A truncated normal continuous random variable.
tukeylambda: A Tukey-Lamdba continuous random variable.
uniform: A uniform continuous random variable.
vonmises: A Von Mises continuous random variable.
vonmises_line: A Von Mises continuous random variable.
wald: A Wald continuous random variable.
weibull_min: Weibull minimum continuous random variable.
weibull_max: Weibull maximum continuous random variable.
wrapcauchy: A wrapped Cauchy continuous random variable.



'''




"""
scipy.stats.probplot returns two tuples of ndarrays:
    First is (osm, osr) or Y and X data:
            Tuple of theoretical quantiles (osm, osr order statistic medians) and ordered responses (osr). osr is simply sorted input x.
    Second is (slope, intercept, r):
            Tuple containing the slope (m) and y-intercept (b) of the straight line and the result of the least-squares fit (r).

for more details visit https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.probplot.html
"""
import pylab 
import scipy.stats as stats #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.probplot.html

vcol = "Nzr"

xy, mbr = stats.probplot(df[f"{vcol}"], dist="norm", fit = True, plot=pylab)
plt.annotate("R^2 = {:.3f}".format(mbr[2]), (xy[0].min(), xy[0].max())) # mbr[2] is the r value, xy[0] is the x-values, i.e. the theoretical quantiles, i.e. line 
plt.title(f"Probability Plot: df[{vcol}]") # Isn't it a Q/Q plot when you are comparing the quantiles?
pylab.show()


fig, ax = plt.subplots(1, 2, figsize=(14, 4))
probplot = sm.ProbPlot(df[f"{vcol}"], dist=lognorm, fit=True)
probplot.ppplot(line='45', ax=ax[0])
probplot.qqplot(line='45', ax=ax[1])
ax[0].set_title(f'P-P Plot: df[{vcol}]')
ax[1].set_title(f'Q-Q Plot: df[{vcol}]')
plt.show()


shape, loc, scale = lognorm.fit(df[f"{vcol}"])
print(shape, loc, scale)
fig, ax = plt.subplots()
stats.probplot(df[f"{vcol}"], fit=True, dist=lognorm, sparams=(shape, loc, scale), plot=ax)
plt.show()


"""Statistics"""

for q in quantiles:
    v = df[vcol].quantile(q)
    print(f"df['{vcol}'] quantile {q} = {v}")

std = df[vcol].std()
var = df[vcol].var()
mean = df[vcol].mean()
cv = std/mean # Coefficient of variation
print(cv)
# =============================================================================
# 
# 
# print('mean=%.3f stdv=%.3f' % (mean(df['z']), std(df['z'])))
# 
# import numpy as np
# import statsmodels.api as sm
# import pylab
# 
# 
# sm.qqplot(df['z'], line='45')
# pylab.show()
# 
# 
# 
# 
# from statsmodels.graphics.gofplots import qqplot
# from matplotlib import pyplot
# 
# qqplot(df['z'], line='s')
# pyplot.show()
# 
# =============================================================================




"""
##############################################################
################ Shapiro-Wilks normality test ################
##############################################################

For N < 5000

In the SciPy implementation of these tests, you can interpret the p value as follows.

p <= alpha: reject H0, not normal.
p > alpha: fail to reject H0, normal.
"""

from scipy.stats import shapiro
# normality test
stat, p = shapiro(df['z'])
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Sample looks Gaussian (fail to reject H0)')
else:
	print('Sample does not look Gaussian (reject H0)')


"""
##############################################################
#################### D’Agostino’s K^2 Test ###################
##############################################################
"""
from scipy.stats import normaltest
# normality test
stat, p = normaltest(df['z'])
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Sample looks Gaussian (fail to reject H0)')
else:
	print('Sample does not look Gaussian (reject H0)')


"""
##############################################################
#################### Anderson-Darling Test ###################
##############################################################

"""
from scipy.stats import anderson

# normality test
result = anderson(df['z'])
print('Statistic: %.3f' % result.statistic)
p = 0
for i in range(len(result.critical_values)):
	sl, cv = result.significance_level[i], result.critical_values[i]
	if result.statistic < result.critical_values[i]:
		print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
	else:
		print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))


"""
##############################################################
#################### Kolgomorov-Smirnov Test ###################
##############################################################

To test for specific distributions.
The null hypothesis is that the data comes from the specified distribution.
If we reject the null, we can conclude that the data do not conform to the tested distribution.

"""
from scipy.stats import kstest
def check_p_val(p_val, alpha):

    if p_val < alpha:
        print('We have evidence to reject the null hypothesis.')
    else:
        print('We do not have evidence to reject the null hypothesis.')

stat, p_val = kstest(df['z'], 'norm')
print('Statistic: \t{:1.2f} \nP-Value: \t{:1.2e}\n'.format(stat, p_val))
check_p_val(p_val, alpha=0.05)


stat, p_val = kstest(df['z'], 'beta', [7, 10])
print('Statistic: \t{:1.2f} \nP-Value: \t{:1.2e}\n'.format(stat, p_val))
check_p_val(p_val, alpha=0.05)







































