# Change_CovFun
change point detection and estimation procedures for the covariance kernel of functional data based on the norms of a generally weighted process of partial sample estimates.

Reference: Change point analysis of covariance functions: a weighted cumulative sum approach, L. Horvath, G. Rice, Y. Zhao (2022) Journal of Multivariate Analysis.

We consider the objective data observations: $X_1(t), \dots, X_N(t)$, $t\in[0,1]$, where each functional observation $X_i(t)$ is a stochastic process with sample path in $L^2([0,1])$. The sequence $X_i(t)$ is assumed to follow the Data Generating Process:

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
\begin{align}\label{model-main}
X_i(t)=\left\{
\begin{array}{ll}
\mu(t)+\eps_i(t),\;\;\;1\leq i \leq k^*
\vspace{.3cm}\\
\mu(t)+\eps_{ i,A}(t), \;\;\;k^*+1\leq i \leq N,
\end{array}
\right.
\end{align}

