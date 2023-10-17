import re
import pandas as pd

# results

results_file = open('./example_results.csv', encoding="UTF-16").readlines()
results = []

for i, lin in enumerate(results_file):
    if "comp" in lin:
        if "WARNING" not in results_file[i+1]:
            molecule = lin.replace('\n','')
            atoms = results_file[i+1].replace('\n','')
            atoms = atoms.split(",")
            atoms = [int(x) for x in atoms]
            atoms2 = results_file[i+2].replace('\n','')
            atoms2 = atoms2.split(",")
            atoms2 = [int(x) for x in atoms2]
        else:
            molecule = lin.replace('\n','')
            atoms = []
            atoms2 = []
        results.append([molecule, atoms, atoms2])

results = pd.DataFrame(results, columns=['molecule','atoms','atoms2'])

# database

database = open('../example/db.smiles').readlines()

for i, lin in enumerate(database):
    database[i] = re.split('\s+',lin)
    database[i] = database[i][0:3]
    database[i][2] = database[i][2].split(",")
    database[i][2] = [int(x) for x in database[i][2] if x != '']

database = pd.DataFrame(database, columns=['molecule','smiles','measured'])

# merge

results = pd.merge(results,database,left_index=True,right_index=True)

corr = 0
semi = 0
incorr = 0
fail = 0

TP = 0
FN = 0

TP2 = 0
FN2 = 0

for i in range(len(results)):

    if results['atoms'][i] == [] and results['atoms2'][i] == []:
        fail = fail + 1
        continue

    atoms = set(results['atoms'][i])
    atoms2 = set(results['atoms2'][i])
    measured = set(results['measured'][i])

    if measured.issubset(atoms): 
        corr = corr+1
    elif measured.issubset(atoms2): 
        semi = semi+1
    else: 
        incorr = incorr+1

    TP += len(atoms.union(measured))
    FN += len(measured - atoms)

    TP2 += len(atoms2.union(measured))
    FN2 += len(measured - atoms2)


print("Number correct = " + str(corr))
print("Number semi-correct = " + str(semi))
print("Number incorrect = " + str(incorr))
print("Number fail = " + str(fail))
print("TPR = " + str(TP/(TP+FN)))
print("TPR2 = " + str(TP2/(TP2+FN2)))

