# -*- coding: utf-8 -*-

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import argparse

parser = argparse.ArgumentParser(description="该脚本用于过滤rdkit无法正确读取的sdf文件")            
parser.add_argument('--filename', '-f', required=True, help='输入旧的sdf文件，格式如：example.sdf')   
parser.add_argument('--outputfile','-o', required=True, help='输入输出sdf文件名称，格式如：output.sdf')
args = parser.parse_args()   
filename = args.filename
outputfile = args.outputfile


suppl = Chem.SDMolSupplier(filename)
writer = Chem.SDWriter(outputfile)
for mol in suppl:
	try:
		AllChem.EmbedMolecule(mol)
		AllChem.MMFFOptimizeMolecule(mol)
		AllChem.AddHs(mol)
		writer.write(mol)
	except Exception as e:
		pass
	continue
writer.close()