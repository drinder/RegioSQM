from rdkit import Chem
from rdkit.Chem import Draw

def draw_sdf(sdf_file_name):
    with Chem.SDMolSupplier('./example/' + sdf_file_name) as suppl:
        ms = [x for x in suppl if x is not None]
    png_file_name = sdf_file_name.replace('.sdf','.png')
    Draw.MolToFile(ms[0],'./images/' + png_file_name)

draw_sdf('comp1+_1-0.sdf')
draw_sdf('comp1+_2-0.sdf')
draw_sdf('comp1+_3-0.sdf')

smiles = 'n1ccc[nH]1'
m = Chem.MolFromSmiles(smiles)
Draw.MolToFile(m,'./images/comp1.png')