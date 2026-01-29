import sys
from rdkit import Chem

print(sys.executable)
print(Chem.MolFromSmiles("CCO"))
