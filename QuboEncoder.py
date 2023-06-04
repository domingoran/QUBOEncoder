import numpy as np


class QEncoder:
    def __init__(
        self,
        Values=[],
        NumberQubits=3,
        EqCoeffs=[[]],
        Lambda=[1],
        Coef=[],
        ScaleFactors=[1],
        BitstringSolution=(),
    ):
        self.NumberVariables = len(EqCoeffs[0])
        self.NumberQubits = NumberQubits
        self.Values = Values
        self.EqCoeffs = EqCoeffs
        self.BitstringSolution = BitstringSolution
        self.Solution = []
        self.EqValues = []
        self.SquareSumConst = 0
        self.Coef = Coef

        if len(Lambda) == 1:
            self.Lambda = Lambda * len(EqCoeffs)
        else:
            self.Lambda = Lambda

        if len(ScaleFactors) != len(EqCoeffs[0]):
            self.ScaleFactors = ScaleFactors * len(EqCoeffs)
        else:
            self.ScaleFactors = ScaleFactors

    def LinearCoef(self, i, l):

        return (self.EqCoeffs[l][i // self.NumberQubits] * self.Coef[i]) * (
            (self.EqCoeffs[l][i // self.NumberQubits] * self.Coef[i])
            - 2 * self.Values[l]
        )

    def CalculateCoef(self):
        if len(self.Coef)!=0:
            self.Coef=[]
        for i in range(self.NumberVariables):
            for j in range(self.NumberQubits):
                self.Coef.append(1 / (2 ** (j + 1) * self.ScaleFactors[i]))

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
                        self.Coef[i]
                        * self.Coef[j]
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
                        self.Coef[i + k * self.NumberQubits]
                        * self.Coef[j + k * self.NumberQubits]
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
        self.CalculateCoef()
        Mat = self.Qconst() + self.Qobj()
        return Mat

    def DecodeSolution(self):
        if len(self.Solution) == 0:
            for i in range(len(self.EqCoeffs[0])):
                s = 0
                for j in range(self.NumberQubits):
                    s += (self.BitstringSolution)[
                        j + i * self.NumberQubits
                    ] * self.Coef[j + i * self.NumberQubits]
                self.Solution.append(s)

    def CalculateSquareSumConstr(self):
        if len(self.Solution) == 0:
            self.DecodeSolution()
        self.SquareSumConst = sum(map(lambda x: x * x, self.Solution))

    def CalculateEquationValues(self):
        if len(self.Solution) == 0:
            self.DecodeSolution()
        for i in range(len(self.EqCoeffs)):
            e = 0
            for j in range(len(self.EqCoeffs[i])):
                e += self.EqCoeffs[i][j] * self.Solution[j]
            self.EqValues.append(round((e - self.Values[i]) ** 2, 5))
