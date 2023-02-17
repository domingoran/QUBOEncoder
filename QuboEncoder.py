import numpy as np


class Qmatrix:
    def __init__(
        self,
        NumberQubits=3,
        ScaleFactors=[1],
        Lambda=[1],
        Values=[],
        EqCoeffs=[[]],
        coef=[],
       
    ):
        
        self.NumberQubits = NumberQubits
        self.Values = Values
        self.EqCoeffs = EqCoeffs
        self.NumberVariables = len(EqCoeffs[0])
        if len(Lambda) != len(EqCoeffs):
            self.Lambda = Lambda * len(EqCoeffs)
        else:
            self.Lambda = Lambda
        self.coef = coef
        if len(ScaleFactors) != len(EqCoeffs[0]):
            self.ScaleFactors = ScaleFactors * len(EqCoeffs)
        else:
            self.ScaleFactors = ScaleFactors

    def LinearCoef(self, i, l):

        return (self.EqCoeffs[l][i // self.NumberQubits] * self.coef[i]) * (
            (self.EqCoeffs[l][i // self.NumberQubits] * self.coef[i])
            - 2 * self.Values[l]
        )

    def Coef(self):

        for i in range(self.NumberVariables):
            for j in range(self.NumberQubits):
                self.coef.append(1 / (2 ** (j + 1) * self.ScaleFactors[i]))

    def Qobj(self):
        Qobj = np.zeros(
            (
                self.NumberVariables * self.NumberQubits,
                self.NumberVariables * self.NumberQubits,
            )
        )

        for l in range(len(self.EqCoeffs)):

            Mat = np.zeros(
                (
                    self.NumberVariables * self.NumberQubits,
                    self.NumberVariables * self.NumberQubits,
                )
            )

            for i in range(self.NumberVariables * self.NumberQubits):
                for j in range(self.NumberVariables * self.NumberQubits):
                    Mat[i, j] = (
                        self.coef[i]
                        * self.coef[j]
                        * self.EqCoeffs[l][i // self.NumberQubits]
                        * self.EqCoeffs[l][j // self.NumberQubits]
                    )

            for i in range(self.NumberVariables * self.NumberQubits):
                for j in range(self.NumberVariables * self.NumberQubits):
                    if i < j:
                        Mat[i, j] = 2 * Mat[i, j]
                    if i > j:
                        Mat[i, j] = 0
                    if i == j:
                        Mat[i, j] = self.LinearCoef(i, l)
            Qobj += self.Lambda[l] * Mat

        return Qobj

    def Qconst(self):

        Mat = np.zeros(
            (
                self.NumberVariables * self.NumberQubits,
                self.NumberVariables * self.NumberQubits,
            )
        )
        subcells = []
        for i in range(self.NumberVariables):
            subcells.append(
                np.zeros(
                    (
                        self.NumberQubits,
                        self.NumberQubits,
                    )
                )
            )

        for k in range(self.NumberVariables):
            for i in range(self.NumberQubits):
                for j in range(self.NumberQubits):
                    subcell = subcells[k]
                    subcell[i, j] = (
                        self.coef[i + k * self.NumberQubits]
                        * self.coef[j + k * self.NumberQubits]
                    )

        for i in range(self.NumberVariables):
            Mat[
                i * self.NumberQubits : (i + 1) * self.NumberQubits,
                i * self.NumberQubits : (i + 1) * self.NumberQubits,
            ] = subcells[i]

        for i in range(self.NumberVariables * self.NumberQubits):
            for j in range(self.NumberVariables * self.NumberQubits):
                if i < j:
                    Mat[i, j] = 2 * Mat[i, j]
                if i > j:
                    Mat[i, j] = 0
        return Mat

    def CalculateMatrix(self):
        if len(self.coef) == 0:
            self.Coef()
        Mat = self.Qconst() + self.Qobj()
        return Mat

