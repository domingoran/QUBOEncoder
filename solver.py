import QuboEncoder as qe
from dwave.system import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite
from minorminer import find_embedding
import dwave.inspector
import pandas as pd
import numpy as np
import datetime
import dimod

n_q = 4
t_an = 20
num_reads = 4000
with open("Results.txt", "a") as f:
    f.write(
        f"----------Run: {datetime.datetime.now()}---n_q = {n_q}---t_an = {t_an}---n_reads = {num_reads}"
    )
    f.write("\n")


arrT = []

# Find values of Opt Forces and save in a list
df = pd.read_excel(r"Dati.xlsx", sheet_name="OptimalForce")
OptForce = list(df.OptForce)

PageNames = [
    "MomentArm_hip_flexion_l",
    "MomentArm_hip_adduction_l",
    "MomentArm_hip_rotation_l",
    "MomentArm_knee_angle_l",
    "MomentArm_ankle_angle_l",
]

# import Data of Momenta, Activations and  ArmsMatrix at time instant i
for i in range(0, 132, 200):
    ArmsMatrix = []
    time = 5.133311 + i * (6.216601 - 5.133311) / 130
    arrT.append(time)
    df = pd.read_excel(r"Dati.xlsx", sheet_name="Momenta")
    Momenta = list(df.iloc[i, 7:12])
    df = pd.read_excel(r"Dati.xlsx", sheet_name="Activations")
    Activations = list(df.iloc[i, 1:11])
    for sheet_name in PageNames:
        df = pd.read_excel(
            r"Dati.xlsx",
            sheet_name=sheet_name,
        )
        Arms = list(df.iloc[i, 1:11])
        ArmsForce = [Arms[i] * OptForce[i] for i in range(len(OptForce))]
        ArmsMatrix.append(ArmsForce)

    Obj = qe.QEncoder(
        NumberQubits=n_q,
        EqCoeffs=ArmsMatrix,
        Values=Momenta,
        Lambda=[abs(1 / (Momentum)) for Momentum in Momenta],
        ScaleFactors=[2, 1.5, 2, 3, 2, 4, 3, 2, 2, 4],
    )
    Mat = Obj.CalculateMatrix()

    solver = DWaveSampler(solver={"topology__type": "pegasus"})

    bqm = dimod.BinaryQuadraticModel.from_qubo(Mat)
    Q = (bqm.to_qubo())[0]
    __, target_edgelist, target_adjacency = solver.structure
    emb = find_embedding(Q, target_edgelist, verbose=1)
    sampler = FixedEmbeddingComposite(solver, emb)
    with open("Embedding.txt", "a") as f:
        f.write(f"----------Run: {datetime.datetime.now()}---n_q = {n_q}")
        f.write(emb)
        f.write(f"Number of physical qubits = {sum(len(chain) for chain in emb.values())} \n")
        f.write("\n")
    results = sampler.sample_qubo(
        Q, num_reads=num_reads, chain_strength=10 * Mat.max(), annealing_time=t_an
    )
    dict_sol = results.first.sample
    keys, values = zip(*dict_sol.items())
    calculationTime=results.info["timing"]

    with open("Temp.txt", "a") as f:
        np.savetxt(f, results)

    Obj.BitstringSolution = values
    Obj.CalculateEquationValues()
    Obj.CalculateSquareSumConstr()

    with open("Results.txt", "a") as f:
        f.write(f"EqValues = {Obj.EqValues} \n")
        f.write(f"Solution = {Obj.Solution} \n")
        f.write(f"Activations = {Activations} \n")
        f.write(f"SquareSumConst = {Obj.SquareSumConst} \n")
        f.write(f"CalculationTime= {calculationTime} \n")


