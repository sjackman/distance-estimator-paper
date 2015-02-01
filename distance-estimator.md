---
title: 'DistanceEstimator: Estimate the distance between two nucleotide sequence fragments using paired-end reads'
author: Shaun D Jackman, Inanc Birol
documentclass: wlpeerj
bibliography: distance-estimator.bib
csl: peerj.csl
header-includes:
  \DeclareMathOperator*{\argmax}{arg\,max}
  \DeclareMathOperator{\MLE}{MLE}
  \keywords{genome sequence assembly, scaffolding, maximum likelihood estimator, distance, gap}
abstract: Paired-end reads can be used to estimate the distance between two nucleotide sequence, called contigs. The difference of the mean of the fragment sizes that span two contigs and the global population of fragment sizes is a trivial but flawed estimate of the size of the gap. Using an empirical fragment size distribution and the maximum likelihood estimator yields more accurate estimates.
---

Methods
================================================================================

To estimate the distance between two contigs, we start by mapping paired-end reads to the two contigs, using software such as ABySS-map, BWA [@li2009fast] or Bowtie2 [@langmead2012fast]. The two contigs are then ordered and oriented to agree with the orientation of the reads. We establish a coordinate system where $-l_1$ and $-1$ are the first and last base of the first contig of length $l_1$, and $0$ and $l_2-1$ are the first and last base of the second contig of length $l_2$.  An observed fragment size, $x_i$, is calculated for each pair by calculating the difference of the mapped position of the first sequenced base of each of the two reads, that is, the outer coordinates of the read pair. This observed fragment size differs from the actual fragment size by the size of the gap between the two contigs, $\theta_0$.

The distribution of fragment sizes of the library is derived empirically by mapping the reads to contigs assembled wihout using paired-end information, or alternatively mapped to a reference genome, and creating a histogram of inferred fragment sizes.

Estimator assuming a normal distribution
------------------------------------------------------------

If we assume that fragment size is normally distributed with mean $\mu$, a
reasonable estimate, $\hat \theta_\text{N}$, of the size of the gap is the
difference between the mean of the population, $\mu$, and the mean of the
sample, $\bar x$.

$$
\hat \theta_\text{N} = \mu - \bar x
$$

Estimator using an empirical distribution
------------------------------------------------------------

The estimate of the distance between the two contigs can be improved by using an empirically derived probability distribution rather than a assuming a normal distribution. Let the probability of observing a fragment of size $x$ selected at random from the population be $f_Z(x)$, and the probability of observing a fragment of size $x$ that spans a gap of size $\theta$ be $f_\theta(x)$. With a sample of $n$ observed fragment sizes, $x_1, \dotsc , x_n$, the likelihood that the two contigs are separated by a distance of $\theta$ bases is $P(X_1 = x_1, \dotsc, X_n = x_n \mid \Theta = \theta)$, or $\mathcal{L}(\theta \mid x_1, \dotsc , x_n)$.

The most likely estimate of the size of the gap between the two contigs is the value $\hat \theta_{\MLE}$ that maximizes the likelihood function, or conveniently, the log likeliehood function, since log is a monotonic transformation.

$$
\begin{aligned}
\hat \theta_{\MLE}
&= \argmax_\theta \mathcal{L}(\theta \mid x_1, \dotsc, x_n)
	= \argmax_\theta \prod_{i=1}^n f_\theta(x_i) \\
&= \argmax_\theta \log \mathcal{L}(\theta \mid x_1, \dotsc, x_n)
	= \argmax_\theta \sum_{i=1}^n \log f_\theta(x_i)
\end{aligned}
$$

Distribution of fragment sizes that span a gap
------------------------------------------------------------

The distribution of observed fragment sizes is the population distribution, $P(Z = x)$, shifted by the size of the gap, $\theta$. Since we can only observe fragments that actually span the gap, we use Bayes' thereom to determine the conditional probability of observing a fragment of size $x$ given that it spans the gap of size $\theta$.

$$
\begin{aligned}
f_{\theta}(x)
&= P(Z = x + \theta \mid \text{fragment spans gap}) \\
&= \frac{ P(\text{fragment spans gap} \mid Z = x + \theta) P(Z = x + \theta) }
{ P(\text{fragment spans gap}) } \\
&\propto P(\text{fragment spans gap} \mid Z = x + \theta) P(Z = x + \theta)
\end{aligned}
$$

Assume the reads are sampled uniformly from the genome between the coordinates $0$ and $g$, where $g$ is the size of the genome, then the position of the left read on the genome $A$ is distributed uniformly with $A \sim U(0, g)$.

$$
P(A = a) = \begin{cases}
\frac{1}{g} & 0 \leq a < g \\
0 & \text{otherwise}
\end{cases}
$$

The probability that a fragment of size $x_i$ spans the gap is the probability that the fragment's left coordinate, $a_i$, falls to the left of the gap, and its right coordinate, $a_i + x_i$, falls to the right of gap. Let $w_\theta(z)$ be the probability that a fragment of observered size $x$ and actual size $z = x + \theta$ spans the gap of size $\theta$.

$$
\begin{aligned}
w_\theta(x + \theta)
&=P(\text{fragment spans gap of size}\ \theta \mid Z = x + \theta) \\
&=P(-l_1 \leq A < 0 \wedge \theta \leq A + Z < l_2 + \theta \mid Z = x + \theta) \\
&= P(-l_1 \leq A < 0 \wedge \theta \leq A + x + \theta < l_2 + \theta) \\
&= P(-l_1 \leq A < 0 \wedge 0 \leq A + x < l_2) \\
&= P(\text{fragment spans gap of size}\ 0 \mid Z = x) \\
&= w_0(x) \\
\end{aligned}
$$

This shows that $w_\theta(x + \theta)$, the probability that a fragment of actual size $x + \theta$ spans the gap, is independent of the size of the gap, $\theta$, and depends only on the observed size of the fragment, $x$. Let $w(x) = w_0(x)$.

$$
\begin{aligned}
w(x) &= w_0(x) \\
&= P(-l_1 \leq A < 0 \wedge 0 \leq A + x < l_2) \\
&= P(-l_1 \leq A < 0 \wedge -x \leq A < l_2 - x) \\
&= P(\max(-l_1, -x) \leq A < \min(0, l_2 - x)) \\
&= \sum_{i=\max(-l_1, -x)}^{\min(0, l_2 - x) - 1} P(A = i) \\
&= \sum_{i=\max(-l_1, -x)}^{\min(0, l_2 - x) - 1} \frac{1}{g} \\
&= \max(0, \min(0, l_2 - x) - \max(-l_1, -x)) \frac{1}{g} \\
&\propto \max(0, \min(0, l_2 - x) - \max(-l_1, -x)) \\
&= \max(0, \min(0, l_2 - x) + \min(l_1, x)) \\
&= \max(0, \min(x, l_1, l_2, l_1 + l_2 - x))
\end{aligned}
$$

Without loss of generality, assume $l_1 \leq l_2$.

$$
\begin{aligned}
w(x) &\propto \begin{cases}
x & 0 \leq x < l_1 \\
l_1 & l_1 \leq x < l_2 \\
l_1 + l_2 - x & l_2 \leq x < l_1 + l_2 \\
0 & \text{otherwise}
\end{cases}
\\
f_\theta(x) &\propto f_Z(x + \theta) w(x)
\\
f_\theta(x) &= \frac{ f_Z(x + \theta) w(x) }
	{ \sum_{j=1}^g f_Z(j + \theta) w(j) }
\end{aligned}
$$

Maximum likelihood estimator
------------------------------------------------------------

We now substitute the distribution of observed fragment sizes, $f_\theta(x)$, into the formula for the maximum likelihood estimator.

$$
\begin{aligned}
\mathcal{L}(\theta \mid x_1, \dotsc, x_n)
&= \prod_{i=1}^n f_\theta(x_i) \\
&= \prod_{i=1}^n \frac{ f_Z(x_i + \theta) w(x_i) }
	{ \sum_{j=1}^g f_Z(j + \theta) w(j) } \\
&= \frac{ \prod_{i=1}^n f_Z(x_i + \theta) w(x_i) }
	{ \left( \sum_{j=1}^g f_Z(j + \theta) w(j) \right) ^n }
\\
\log \mathcal{L}(\theta \mid x_1, \dotsc, x_n)
&= \sum_{i=1}^n \log f_Z(x_i + \theta)
	+ \sum_{i=1}^n \log w(x_i)
	- n \log \sum_{j=1}^g f_Z(j + \theta) w(j)
\\
\hat \theta_{\MLE}
&= \argmax_\theta \log \mathcal{L}(\theta \mid x_1, \dotsc, x_n) \\
&= \argmax_\theta \left[ \sum_{i=1}^n \log f_Z(x_i + \theta)
	+ \sum_{i=1}^n \log w(x_i)
	- n \log \sum_{j=1}^g f_Z(j + \theta) w(j) \right] \\
&= \argmax_\theta \left[ \sum_{i=1}^n \log f_Z(x_i + \theta)
	- n \log \sum_{j=1}^g f_Z(j + \theta) w(j) \right]
\end{aligned}
$$

Finding the value of $\theta$ that maximizes the likelihood function is an optimization problem. When the fragment size of the sequencing library is small, the range of possible values of $\theta$ is small, and it is reasonable to calculate every value of $\mathcal{L}(\theta)$ to find the maximum.

Conclusion
================================================================================

This distance estimation algorithm is implemented by the ABySS assembly software [@simpson2009abyss] in the utility DistanceEst, which takes as input the distribution of fragment sizes of the sequencing library and a SAM file [@li2009sequence] of paired-end reads that map to different contigs.

Acknowledgements
================================================================================

The distance estimator described here is derived from the estimator implemented by Jared Simpson in the first release of the software ABySS.

References
================================================================================
