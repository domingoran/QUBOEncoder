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
import os
import pathlib

def scalingCalculator(lastScaling,max,n_q):
    return  1/(1/(lastScaling*2**(n_q))+max)

n_q = 4
t_an = 20
num_reads = 4000
L=0.1
n_sim=0
n_scaling=2

OptForce = [1313.18360656, 557.1147541 ,2191.74098361,5148.76721311,2747.82295082]

for n_s in range(n_scaling+1):
    solutions={}
    path="Results/L"+str(L)+"/sim"+str(n_sim)
    pathDir=path+'/nScal'+str(n_s)
    pathlib.Path(pathDir).mkdir(parents=True, exist_ok=True) 

    if (os.path.exists(path+"/scaling.json") == False):
        scaling={}
        scaling[0]= [1.0, 1.0, 1.0, 1.0, 1.0]
        with open(path+'/scaling.json', 'w') as f:
            json.dump(scaling, f)

    g = open(path+'/scaling.json')
    scaling=json.load(g)
    sf=scaling[str(n_s)]
    for i in range(0,150,100):
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
            ScaleFactors=sf, 
            Coef=[]
        )
    
        Mat = Obj.CalculateMatrix()
        
        solver = DWaveSampler(solver={"topology__type": "pegasus"})
        

        bqm = dimod.BinaryQuadraticModel.from_qubo(Mat)
        Q = (bqm.to_qubo())[0]
        __, target_edgelist, target_adjacency = solver.structure
        emb = find_embedding(Q, target_edgelist, verbose=1)
        sampler = FixedEmbeddingComposite(solver, emb)

        filename = pathDir+"/Result_lambda_"+str(L)+"_sim_"+str(n_sim)+"_nq_"+str(n_q)

        results = sampler.sample_qubo(
            Q, num_reads=num_reads, chain_strength=10 * Mat.max(), annealing_time=t_an
        )
        dict_sol = results.first.sample
        keys, values = zip(*dict_sol.items())
        
        df1=pd.DataFrame(list(map(np.ravel, results.record["sample"][0:100])))
        df2 = pd.DataFrame({'Energy': results.record["energy"][0:100], 'num_occurences': results.record["num_occurrences"][0:100]})
        df=pd.concat([df1, df2], axis="columns")
        df.to_csv(pathDir+"/Energies_lambda_"+str(L)+"_sim_"+str(n_sim)+"_nq_"+str(n_q)+'_t_'+str(i)+'.csv')

        Obj.BitstringSolution = values
        Obj.CalculateEquationValues()
        Obj.CalculateSquareSumConstr()

        solutions[i/10]=Obj.Solution

        n_p_q = []
        for l in range(len(emb)):   
            n_p_q.append(len(emb[l]))

        result = ({'t': i, 
                'constraint': Obj.EqValues, 
                'sol': Obj.Solution, 
                'act': Activations, 
                'fun': Obj.SquareSumConst, 
                'emb': n_p_q,
                'phys_qubits': sum(len(chain) for chain in emb.values())
                } | results.info["timing"])

        pd.DataFrame([result]).to_csv(filename+".csv", mode= 'a', header = not os.path.exists(filename+".csv"))

    scalingArr=[]
    for n_act in range(5):
        maxSol=0
        for el in solutions.values():
            if el[n_act]>maxSol:
                maxSol=el[n_act]
        scalingArr.append(scalingCalculator(sf[n_act],maxSol,n_q))

            
    scaling.update({str(n_s+1):scalingArr})
    with open(path+'/scaling.json', 'w') as f:
        json.dump(scaling, f)
