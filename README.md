---
title: "QUBO Matrix encoder"
author: "Domingo Ranieri"
V0 date: "15 February 2023"
V1 date: "04 April 2023"
---

# QUBO Matrix encoder

QUBO matrix is a useful format to submit problems on a D-Wave quantum computer.
This modulo encodes a specific problem into the QUBO format.
If we have a system of equations with more variables than equations:

$$
\begin{cases} 
-4x_1+2x_2+5x_3-2x_4=9\\\\
7x_1+3x_2-2x_3+1x_4=4
\end{cases}
$$

This system has an infinite set of solutions.\
Imposing these constraints on the solution:

$$
\min\sum_{i=1}^4(x_i^2)
$$

$$
0 \leq x_i \leq 1
$$

it becomes an optimization problem with a unique solution.
In general, the final matrix is:

$$
Q=Q_{const}+ \sum_i\lambda_i Q_{(obj)i}
$$

where:

* $Q_{const}$ is the matrix that encodes the first constraint
* $\lambda_i$ are Lagrangian multipliers that must be tuned
* $Q_{(obj)i}$ are the matrices that encode each equation

The following code can be used to encode the problem assuming that $\lambda_1=0.5$ and $\lambda_2=0.3$ :

```python
import QuboEncoder as qe
Q=qe.QEncoder(EqCoeffs=[[-4,2,5,-2],[7,3,-2,1]], Values=[9,4], Lambda=[0.5,0.3])
QMatrix=Q.CalculateMatrix()
```
If $\lambda_i=\lambda=0.5 \forall i$, the code is modified setting Lambda=[0.5].\
The default value is Lambda=[1]

To satisfy the second constrain, each variable is encoded using a binary sum: 

$$
x_i=\frac{1}{s}\sum_{j=1}^{n_q}\frac{q_j}{2^j}\\
$$

* $n_q$ is the number of qubits used and this is the same $\forall$ $x_i$. By default $n_q=3$.
* $s$ is a scale factor that can be set to reduce the range of values of $x_i$. For each variable can be set a different value of $s$. By defoult each $s=1$ 

From the previous example, to encode each variable using 4 qubits and set the scale factors respectively equal to 2,1,4 and 10 the code becomes:


```python
import QuboEncoder as qe
Q=qe.QEncoder(EqCoeffs=[[-4,2,5,-2],[7,3,-2,1]], Values=[9,4], Lambda=[0.5,0.3], NumberQubits=4, ScaleFactors=[2,1,4,10])
QMatrix=.CalculateMatrix()
```
## Postprocessing results


When the quantum computation is finished, the bitstring solution can be stored and processed with some methods.

If `results` contains the output of the quantum annealing, the bitstring solution can be extracted and stored in this way:

```python
dict_sol = results.first.sample
keys, values = zip(*dict_sol.items())
Q.BitstringSolution=values
```

The method `DecodeSolution`. can be used to pass from the bitstring solution to the $x_i$ values $\in$ [0,1]. The result is stored in the class attribute `Solution`.

It is possible to calculate the value of the first imposed constraint using the method `CalulateSquareSumConstr` and the result is saved in the class attribute `SquareSumConst`.

The method `CalulateEquationValues` calculate the square between the $L.H.S$ and the $R.H.S.$ for each equation in the system. This value should be 0 for the correct solution of the system. The values are stored in the class attribute list `EqValues`.
