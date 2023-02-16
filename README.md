---
title: "QUBO Matrix encoder"
author: "Domingo Ranieri"
date: "15 February 2023"
output: html_document
---

## QUBO Matrix encoder
---
QUBO matrix is a useful format to submit problems on D-Wave quantum computer.\
This modulo encode a specific problem into the QUBO format.\
If we have a system of equations with more variables than equations:

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r equation}
eq <-  noquote(paste(expression("2x+5y-2z+2p=9\\\\3x-2y+1z-3p=34\\\\-3x+3y+2z+4p=33\\\\2x+3y+4z+5p=125")))
```

$$
\begin{cases} `r eq` \end{cases}
$$

$$
\begin{cases} 
-4x_1+2x_2+5x_3-2x_4=9\\\\
7x_1+3x_2-2x_3+1x_4=4
\end{cases}
$$
we know that we have infinite set of solutions.\
Imposing these constraints on solution:
$$
\min\sum_{i=1}^4(x_i^2)\\
0 \leq x_i \leq 1
$$

it becames an optimization problem with a unique solution.\
The final matrix in general is:
$$
Q=Q_{const}+ \sum_i\lambda_i Q_{(obj)i}
$$
where:
* $Q_{const}$ is the matrix that encode the first constraint
* $\lambda_i$ are Lagrangian multipliers that must be tuned
* $Q_{(obj)i}$ are the matrices that encode each equation

The following code can be used to encode the problem assuming that $\lambda_1=0.5$ and $\lambda_2=0.3$ :

```python
import QuboEncoder as qe
Q=qe.Qmatrix(EqCoeffs=[[-4,2,5,-2],[7,3,-2,1]],  Values=[9,4],Lambda=[0.5,0.3]).CalculateMatrix()
```
If $\lambda_i=\lambda=0.5 \forall i$, the code is modified setting Lambda=[0.5].\
The default value is Lambda=[1]

To satisfy the second constrain, each variable is encoded using a binary sum: 

$$
x_i=\frac{1}{s}\sum_{j=1}^{n_q}\frac{q_j}{2^j}\\
$$
* $n_q$ is the number of qubits used and this is the same $\forall$ $x_i$. By defoult $n_q=3$.
* $s$ is a scale factor that can be set to reduce the range of values of $x_i$. For each variables can be set a different value of $s$. By defoult each $s=1$ 

From the previous example, to encode each variables using 4 qubits and set the scale factors respectively equal to 2,1,4 and 10 the code becames:


```python
import QuboEncoder as qe
Q=qe.Qmatrix(EqCoeffs=[[-4,2,5,-2],[7,3,-2,1]],  Values=[9,4],Lambda=[0.5,0.3],NumberQubits=4, scaleFactors=[2,1,4,10]).CalculateMatrix()
```

<!-- If $\lambda_i=\lambda, \forall i$, we get:
$$
Q_{obj}=\sum_i Q_{(obj)i}
$$
and 
$$
Q=Q_{const}+ \lambda Q_{obj}
$$ -->
