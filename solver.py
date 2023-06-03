import QuboEncoder as qe
from dwave.system import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite
from minorminer import find_embedding
import dwave.inspector
import pandas as pd
import numpy as np
import datetime
import dimod
import json


solutions={}


n_q = 4
t_an = 20
num_reads = 4000
L=10
with open("Result10.txt", "a") as f:
    f.write(
        f"----------Run: {datetime.datetime.now()}---n_q = {n_q}---t_an = {t_an}---n_reads = {num_reads}---L={L}"
    )
    
OptForce = [1313.18360656, 557.1147541 ,2191.74098361,5148.76721311,2747.82295082]

# import Data of Momenta, Activations and  ArmsMatrix at time instant i
for i in range(0,150,190):
    ArmsMatrix = []
    df = pd.read_csv(r"SystemSO_KE.csv")
    Momenta = list(df.iloc[i, 5:6])
    Arms=list(df.iloc[i, :5])
    Activations = list(df.iloc[i, 6:-1])
    ArmsForce = [Arms[i] * OptForce[i] for i in range(len(OptForce))]
    ArmsMatrix.append(ArmsForce)

    Obj = qe.QEncoder(
        NumberQubits=n_q,
        EqCoeffs=ArmsMatrix,
        Values=[Momenta[i] -150 for i in range(len(Momenta))],
        # Lambda=[10*abs(1 / (Momenta[0]))],
        Lambda=[L],
        ScaleFactors=[1,1,1,1,1], #[40,40,1.5,1,1.5],
    )
    Mat = Obj.CalculateMatrix()
    solver = DWaveSampler(solver={"topology__type": "pegasus"})
    

    bqm = dimod.BinaryQuadraticModel.from_qubo(Mat)
    Q = (bqm.to_qubo())[0]
    __, target_edgelist, target_adjacency = solver.structure
    emb = find_embedding(Q, target_edgelist, verbose=1)
    sampler = FixedEmbeddingComposite(solver, emb)
    with open("Embedding10.txt", "a") as f:
        f.write(f"----------Run: {datetime.datetime.now()}---n_q = {n_q} \n")
        f.write(str(emb))
        f.write(f"Number of physical qubits = {sum(len(chain) for chain in emb.values())} \n")
        f.write("\n")
    results = sampler.sample_qubo(
        Q, num_reads=num_reads, chain_strength=10 * Mat.max(), annealing_time=t_an
    )
    dict_sol = results.first.sample
    keys, values = zip(*dict_sol.items())
    calculationTime=results.info["timing"]

    df1=pd.DataFrame(list(map(np.ravel, results.record["sample"][0:100])))

    df2 = pd.DataFrame({'Energy': results.record["energy"][0:100], 'num_occurences': results.record["num_occurrences"][0:100]})
    df=pd.concat([df1, df2], axis="columns")
    df.to_csv('allResults/L10/'+str(datetime.datetime.now())+"-"+str(n_q)+"-"+str(t_an)+"-"+str(L)+"-"+str(i)+'.csv')

    Obj.BitstringSolution = values
    Obj.CalculateEquationValues()
    Obj.CalculateSquareSumConstr()

    solutions[i/10]=Obj.Solution

    with open("Result10.txt", "a") as f:
        f.write("\n")
        f.write(f"Time = {i} \n")
        f.write(f"EqValues = {Obj.EqValues} \n")
        f.write(f"Solution = {Obj.Solution} \n")
        f.write(f"Activations = {Activations} \n")
        f.write(f"SquareSumConst = {Obj.SquareSumConst} \n")
        f.write(f"CalculationTime= {calculationTime} \n")

with open('solutions.json', 'w') as f:
    json.dump(solutions, f)

