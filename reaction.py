# -*- coding: utf-8 -*-

#############################
#Load module
#############################
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import argparse
import sys

#############################
#Parsing parameters
#############################
parser = argparse.ArgumentParser(description="该脚本用于生成基于化学反应的组合化合物库构建  作者：grosetta  微信公众号：grosetta")            
parser.add_argument('--filename', '-f', required=True, help='输入用于“反应”的底物化合物库，格式如：example.sdf')   
parser.add_argument('--outputfile','-o', required=True, help='定义输出的组合库名称，格式如：output.sdf')
parser.add_argument('--mode','-m', required=True, type=int, help='化学反应的类型（输入反应类型对应的数字）： 1. 生成酰胺反应 2. 酯化反应 3. 狄尔斯-阿尔德反应')           
args = parser.parse_args()   
filename = args.filename
outputfile = args.outputfile
mode = args.mode

#############################
#read molcules
#############################
suppl = Chem.SDMolSupplier(filename)

def print_list(lst,lis=[]):
	for x in lst:
		if type(x) is list:
			print_list(x)
		else:
			lis.append(x)
	return lis

def print_tuple(lst,lis=[]):
		
		for x in lst:
			if type(x) is tuple:
				lis.append(x[0])
		return lis

def reaction(Substrate1, Substrate2, Equation):
	#############################
	#match pattern1
	#############################
	patt1 = Chem.MolFromSmarts(Substrate1)
	matches_patt1 = []
	for mol in suppl:
		if mol.HasSubstructMatch(patt1,useChirality=True):   
			matches_patt1.append(mol)
	
	#############################
	#match pattern2
	#############################
	patt2 = Chem.MolFromSmarts(Substrate2)
	matches_patt2 = []
	for mol in suppl:
		if mol.HasSubstructMatch(patt2,useChirality=True):   
			matches_patt2.append(mol)
	num1 = len(matches_patt1)
	num2 = len(matches_patt2)
	end1 = num1
	end2 = num2
	if num1 > 50:
		print("匹配到的用于反应的底物1的数量过多，总数为：")
		print(num1)
		judge = input("你想要自定义参与反应的底物1的数量吗？输入Y表示确认，输入N默认用全部的底物，这可能会耗费很长时间 (Y or N)?  ")

	if judge == "Y" :
		end1 = int(input("请输入用于反应的底物1的数量:"))
	elif judge == "N":
		print("")
	else:
		print("错误！！！请输入Y或N！！！")
		sys.exit()

	if end1 > num1:
		print('错误！超过底物总数！')
		sys.exit()
		
	if num2 > 50:
		print("匹配到的用于反应的底物2的数量过多，总数为：")
		print(num2)
		judge = input("你想要自定义参与反应的底物2的数量吗？输入Y表示确认，输入N默认用全部的底物，这可能会耗费很长时间 (Y or N)?  ")

	if judge == "Y":
		end2 = int(input("请输入用于反应的底物2的数量:"))
	elif judge == "N":
		print("")
	else:
		print("错误！！！请输入Y或N！！！")
		sys.exit()

	if end2 > num2:
		print('错误！超过底物总数！')
		sys.exit()
	
	#############################
	#Define the reaction equation
	#############################
	rxn = AllChem.ReactionFromSmarts(Equation)
	result = [ ]
	for i in matches_patt1[0:end1]:
		for j in matches_patt2[0:end2]:
			try:
				ps = rxn.RunReactants((i,j))
				ps_list = list(ps)
				result.append(ps_list)
			except Exception as e:
				pass
			continue
	out = print_list(result)
	final_out = print_tuple(out)
	return final_out
	
	
if mode == 1:
	final_out = reaction("[CX3](=O)[OX2H1]", "[NX3;H2,H1;!$(NC=O)]", "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
if mode == 2:
	final_out = reaction("[CX3](=O)[OX2H1]", "[O;H1,H0;!$(OC=O)]", "[C:1](=[O:2])-[OD1].[O!H0:3]>>[C:1](=[O:2])[O:3]")
if mode == 3:
	final_out = reaction('[#8]-[#6]=[#6]', '[#6]=[#6]-[#6](-[#7])=[#6]', '[C:1]=[C:2].[C:3]=[*:4][*:5]=[C:6]>>[C:1]1[C:2][C:3][*:4]=[*:5][C:6]1')

#############################
#Output result
#############################
writer = Chem.SDWriter(outputfile)
for mol in final_out:
	try:
		AllChem.EmbedMolecule(mol, useRandomCoords=True)
		AllChem.MMFFOptimizeMolecule(mol)
		AllChem.AddHs(mol)
		writer.write(mol)
	except Exception as e:
		pass
	continue
writer.close()