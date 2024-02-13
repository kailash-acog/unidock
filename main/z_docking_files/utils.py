import os
import re
import csv 
import sys
import time
import math
import torch
import random
import shutil
import psutil
import string
import logging
import subprocess
import pandas as pd
import concurrent.futures
import ipywidgets as widgets
import multiprocessing as mp

from glob import glob
from typing import Optional, List
from IPython.display import display
from rich.progress import Progress
from openbabel import openbabel, pybel
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw 
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool, cpu_count

try:
    from protonator import protonator
    has_protonator = True
except ImportError:
    protonator, has_protonator = None, False




##############################################################################################################################
""" Check GPU Availability """

def check_availability():
    if "CUDA_VISIBLE_DEVICES" not in os.environ:
        os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    if torch.cuda.is_available():
        device = torch.device("cuda")
        gpu_info = os.popen('nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits').readlines()
        gpu_available = 100 - int(gpu_info[0].strip())
        gpu_result = f"\033[1m\033[34mGPU availability in own5: \033[91m{gpu_available:.2f}%\033[0m"
    else:
        device = torch.device("cpu")
        gpu_result = 'GPU is not available, using CPU instead'

    cpu_percentage = psutil.cpu_percent()
    cpu_available = 100 - cpu_percentage
    cpu_result = f"\033[1m\033[34mCPU availability in own5: \033[91m{cpu_available:.2f}%\033[0m"
    
    print(gpu_result)
    print(cpu_result)
    return device



##############################################################################################################################
""" Convert SMILES to SDF """

def read_smi_file(filename: str, i_from: int, i_to: int) -> List[Chem.Mol]:
    mol_list = []
    with open(filename, 'r') as smiles_file:
        for i, line in enumerate(smiles_file):
            if i_from <= i < i_to:
                tokens = line.split()
                smiles = tokens[0]
                mol_list.append(Chem.MolFromSmiles(smiles))
    return mol_list


def get_structure(mol: Chem.Mol, num_conformations: int, index: int) -> Optional[Chem.Mol]:
    try:
        if has_protonator:
            mol = protonator(mol)

        mol = Chem.AddHs(mol)
        new_mol = Chem.Mol(mol)

        conformer_energies = []
        AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        conformer_energies = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000, nonBondedThresh=100.0)

        if index == 0:
            i = conformer_energies.index(min(conformer_energies))
        elif index > 0:
            i = index - 1
        else:
            raise ValueError("index cannot be less than zero.")

        new_mol.AddConformer(mol.GetConformer(i))
        return new_mol
    except ValueError as e:
        print(f"Error processing molecule: {e}")
        return None

def molecules_to_structure(population: List[Chem.Mol], num_conformations: int, index: int, num_cpus: int):
    with mp.Pool(num_cpus) as pool:
        args = [(p, num_conformations, index) for p in population]
        generated_molecules = pool.starmap(get_structure, args)

        names = [''.join(random.choices(string.ascii_uppercase + string.digits, k=6)) for _ in generated_molecules]
        return generated_molecules, names


def molecule_to_sdf(mol: Chem.Mol, output_filename: str, name: Optional[str] = None):
    if name is not None:
        mol.SetProp("_Name", name)
    writer = Chem.SDWriter(output_filename)
    writer.write(mol)
    writer.close()


def process_row(row, output_sdf, num_conformations, idx_conformer):
    smiles, mol_name = row['SMILES'], row['Name']
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mol = get_structure(mol, num_conformations, idx_conformer)
        if mol is not None:
            sdf_filename = os.path.join(output_sdf, f"{mol_name}.sdf")
            molecule_to_sdf(mol, sdf_filename, name=mol_name)


def convert_smiles_to_sdf_parallel(df, output_sdf, num_conformations, idx_conformer=0):
    total = len(df)
    max_workers = os.cpu_count()  
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_row, row, output_sdf, num_conformations, idx_conformer) for _, row in df.iterrows()]

        for future in concurrent.futures.as_completed(futures):
            future.result()

    print(f"\033[1m\033[34mSMILES to SDF conversion completed and files saved in folder: \033[91m{output_sdf}\033[0m")




##############################################################################################################################
""" Pass Correct PDBQT files for Docking """

def check_pdbqt_files(folder_name, input_smiles):

    # sys.stderr = open(os.devnull, 'w')
    # sys.stdout = open(os.devnull, 'w')

    logging.getLogger("rdkit").setLevel(logging.ERROR)
    
    def process_smiles_files(input_smiles_files):
        output_file = os.path.join(folder_name, "pipeline_files/smiles.txt")
        if os.path.exists(output_file):
            os.remove(output_file)

        with open(output_file, 'a') as output:
            for filename in os.listdir(input_smiles_files):
                input_file_path = os.path.join(input_smiles_files, filename)
                if os.path.isfile(input_file_path):
                    with open(input_file_path, 'r') as input_file:
                        file_content = input_file.read().strip()
                        output.write(f'{file_content}\n')

        df = pd.read_csv(output_file, header=None, delimiter='\t', names=['obabel_SMILES', 'Name'])
        df = df[['Name', 'obabel_SMILES']]
        return df

    df1 = pd.read_csv(input_smiles)
    df1 = df1.sort_values(by='Name').reset_index(drop=True)

    input_smiles_files = os.path.join(folder_name, "pipeline_files/4_smiles")
    df2 = process_smiles_files(input_smiles_files)
    df2 = df2.sort_values(by='Name').reset_index(drop=True)

    df3 = pd.merge(df1, df2, on='Name', how='inner')

    def process_df(folder_name, df3):
        df4 = pd.DataFrame({'Name': df3['Name']})
        df4['similarity_score'] = 0.0  

        def calculate_similarity(row):
            m1 = Chem.MolFromSmiles(row['SMILES'])
            m2 = Chem.MolFromSmiles(row['obabel_SMILES'])

            if m1 is not None and m2 is not None:
                invgen = AllChem.GetMorganFeatureAtomInvGen()
                ffpgen = AllChem.GetMorganGenerator(radius=2, atomInvariantsGenerator=invgen)

                ffp1 = ffpgen.GetSparseCountFingerprint(m1)
                ffp2 = ffpgen.GetSparseCountFingerprint(m2)

                similarity_score = DataStructs.DiceSimilarity(ffp1, ffp2)
                return similarity_score
            else:
                return None

        df4['similarity_score'] = df3.apply(calculate_similarity, axis=1)
        df4 = df4[df4['similarity_score'] == 1].drop(['similarity_score'], axis=1)

        file_path = f"{folder_name}/pipeline_files/1_compounds_for_docking.csv"
        df4.to_csv(file_path, index=False)
        return df4

    df4 = process_df(folder_name, df3)


def copy_correct_pdbqt_files(folder_name, input_smiles):
    
    all_pdbqt_files = os.path.join(folder_name, "pipeline_files/3_pdbqt")
    compounds_to_be_dock = os.path.join(folder_name, "pipeline_files/1_compounds_for_docking.csv")
    output_dir = os.path.join(folder_name, "pipeline_files/5_pdbqt_for_docking")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df1 = pd.read_csv(input_smiles)
    df2 = pd.read_csv(compounds_to_be_dock)
    filtered_out = len(df1) - len(df2)
    
    compounds_df = pd.read_csv(compounds_to_be_dock)

    for compound_name in compounds_df['Name']:
        input_file_path = os.path.join(all_pdbqt_files, f"{compound_name}.pdbqt")
        output_file_path = os.path.join(output_dir, f"{compound_name}.pdbqt")

        if os.path.exists(input_file_path):
            shutil.copy(input_file_path, output_file_path)

    print(f"\033[1m\033[34mCompounds filtered out using Dice Similarity: \033[91m{filtered_out}\033[0m")




##############################################################################################################################
""" Create a ligands paths text file """

def create_ligands_path_batchwise(folder_name, batch_size=20000):
    output_pdbqt = os.path.join(folder_name, "pipeline_files/5_pdbqt_for_docking")

    def chunk_list(input_list, chunk_size):
        return [input_list[i:i + chunk_size] for i in range(0, len(input_list), chunk_size)]

    def get_pdbqt_files(input_path):
        return [file for file in os.listdir(input_path) if file.endswith(".pdbqt")]

    pdbqt_files = get_pdbqt_files(output_pdbqt)

    ligand_batches = chunk_list(pdbqt_files, batch_size)

    for i, ligand_batch in enumerate(ligand_batches):
        batch_ligands_path = os.path.join(output_pdbqt, "..", f"unidock_pdbqt_batch_{i+1}.txt")
        with open(batch_ligands_path, "w") as batch_file:
            batch_file.write('\n'.join(os.path.join(output_pdbqt, file) for file in ligand_batch))

    batch_files = [f"ligands_batch_{i}.txt" for i in range(len(ligand_batches))]
    return len(batch_files)




##############################################################################################################################
""" Extract Affinity Values """

def affinity_from_pdbqt_files(folder_name):
    ligands_pdbqt_out = os.path.join(folder_name, "pipeline_files/6_pdbqt_out")
    results = []
    for filename in os.listdir(ligands_pdbqt_out):
        file_path = os.path.join(ligands_pdbqt_out, filename)

        if os.path.isfile(file_path) and filename.endswith(".pdbqt"):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                affinity_line = lines[1]
                affinity_value = float(affinity_line.split()[3])
                name = filename.replace('_out.pdbqt', '')
                results.append({'Name': name, 'Affinity': affinity_value})

    output_file = os.path.join(folder_name, 'pipeline_files/2_extract_affinity_from_pdbqt.csv')
    with open(output_file, 'w', newline='') as csv_file:
        fieldnames = ['Name', 'Affinity']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\033[1m\033[34mAffnity values extracted and saved in folder: \033[91m{output_file}\033[0m")




#############################################################################################################################
""" Extract Compounds Based on Affinity threshold """

def extraction_based_on_threshold(folder_name, threshold):
    source_dir = os.path.join(folder_name, "pipeline_files/6_pdbqt_out")
    affinity_score_path = os.path.join(folder_name, "pipeline_files/2_extract_affinity_from_pdbqt.csv")

    df = pd.read_csv(affinity_score_path)

    destination_dir = os.path.join(folder_name, "pipeline_files/7_pdbqt_out_threshold")
    os.makedirs(destination_dir, exist_ok=True)

    output_file_path = os.path.join(folder_name, "pipeline_files/3_compounds_for_posebusters.csv")

    if threshold == 'dynamic':
        dynamic_threshold = df['Affinity'].mean() - df['Affinity'].std()
        dynamic = df[df['Affinity'] < dynamic_threshold]
        dynamic.to_csv(output_file_path, index=False)

        for index, row in dynamic.iterrows():
            compound_name = row['Name'].split('_out')[0]
            source_file = os.path.join(source_dir, f"{compound_name}_out.pdbqt")
            destination_file = os.path.join(destination_dir, f"{compound_name}_out.pdbqt")
            shutil.copy(source_file, destination_file)

    elif isinstance(threshold, (float, int)):
        static = df[df['Affinity'] < threshold]
        static.to_csv(output_file_path, index=False)

        for index, row in static.iterrows():
            compound_name = row['Name'].split('_out')[0]
            source_file = os.path.join(source_dir, f"{compound_name}_out.pdbqt")
            destination_file = os.path.join(destination_dir, f"{compound_name}_out.pdbqt")
            shutil.copy(source_file, destination_file)

    print("\033[1m\033[34mCompounds Extracted based on threshold value\033[0m".format(output_file_path))



#############################################################################################################################
""" Extracted Model 1 content """

def copy_content(file_path, output_folder):
    with open(file_path, 'r') as file:
        content = file.read()
        endmdl_index = content.find("ENDMDL")
        endmdl_content = content[:endmdl_index + len("ENDMDL")]
        output_file_path = os.path.join(output_folder, os.path.basename(file_path))
        with open(output_file_path, 'w') as output_file:
            output_file.write(endmdl_content)
        

def extract_model1(folder_name):
    input_folder = os.path.join(folder_name, "pipeline_files/7_pdbqt_out_threshold")
    output_folder = os.path.join(folder_name, "pipeline_files/8_pdbqt_out_threshold_m1")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file_name)
        if os.path.isfile(file_path) and file_name.endswith(".pdbqt"):
            copy_content(file_path, output_folder)

    print(f"\033[1m\033[34mExtracted Model_1 content and saved in folder: \033[91m{output_folder}\033[0m")



##############################################################################################################################
""" Process PoseBusters Output file """

def process_pb_csv(folder_name):
    pb_result = os.path.join(folder_name, 'pipeline_files', '5_pb_out.csv')
    pb = pd.read_csv(pb_result)
    pb = pb.drop('file', axis=1)
    pb = pb.drop(pb.index[1::2])
    pb = pb.rename(columns={'molecule': 'Name'})
    pb['passes'] = pb.iloc[:, 1:].eq('True').sum(axis=1)
    pb.to_csv(os.path.join(folder_name, 'pipeline_files', '6_pb_out.csv'), index=False)


def final_output(folder_name, input_smiles, passes):
    posebusters_path = os.path.join(folder_name, 'pipeline_files', '6_pb_out.csv')
    df = pd.read_csv(posebusters_path)
    
    df1 = df[df['passes'] >= passes]

    affinity_path = os.path.join(folder_name, 'pipeline_files', '2_extract_affinity_from_pdbqt.csv')
    df2 = pd.read_csv(affinity_path)
    df2 = df2.sort_values(by='Affinity').reset_index(drop=True)

    df_temp = pd.merge(df1, df2, on='Name', how='left')

    df3 = pd.read_csv(input_smiles)

    df4 = pd.merge(df_temp, df3, on='Name', how='left')
    df4 = df4[['Name', 'SMILES', 'Affinity']]
    df4 = df4.sort_values(by='Affinity').reset_index(drop=True)
    output_path = os.path.join(folder_name, 'output.csv')
    df4.to_csv(output_path, index=False)

    filtered_count = len(df) - len(df1)
    print(f"\033[1m\033[34mCompounds filtered out by PoseBusters: \033[91m{filtered_count}\033[0m")



##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


def copy_pdbqt_files_flex(df, folder_name, flex_folder_name):
    source_dir = os.path.join(folder_name, "pipeline_files/5_pdbqt_for_docking")
    destination_dir = os.path.join(flex_folder_name, "pipeline_files/1_pdbqt_in_flex")
    os.makedirs(destination_dir, exist_ok=True)

    for name in df['Name']:
        file_name = name + ".pdbqt"
        source_path = os.path.join(source_dir, file_name)
        destination_path = os.path.join(destination_dir, file_name)
        shutil.copy(source_path, destination_path)


def create_ligands_path_flex(folder_name):
    output_pdbqt = os.path.join(folder_name, "pipeline_files/1_pdbqt_in_flex")

    def get_pdbqt_files(input_path):
        return [file for file in os.listdir(input_path) if file.endswith(".pdbqt")]

    pdbqt_files = get_pdbqt_files(output_pdbqt)

    ligands_path = os.path.join(output_pdbqt, "..", "unidock_pdbqt_path.txt")
    with open(ligands_path, "w") as batch_file:
        batch_file.write('\n'.join(os.path.join(output_pdbqt, file) for file in pdbqt_files))


def affinity_from_pdbqt_files_flex(folder_name):
    ligands_pdbqt_out = os.path.join(folder_name, "pipeline_files/2_pdbqt_out_flex")
    results = []
    for filename in os.listdir(ligands_pdbqt_out):
        file_path = os.path.join(ligands_pdbqt_out, filename)

        if os.path.isfile(file_path) and filename.endswith(".pdbqt"):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                affinity_line = lines[1]
                affinity_value = float(affinity_line.split()[3])
                name = filename.replace('_out.pdbqt', '')
                results.append({'Name': name, 'Affinity': affinity_value})

    output_file = os.path.join(folder_name, 'pipeline_files/1_extract_affinity_from_pdbqt.csv')
    with open(output_file, 'w', newline='') as csv_file:
        fieldnames = ['Name', 'Affinity']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\033[1m\033[34mAffnity values extracted and saved in folder: \033[91m{output_file}\033[0m")


def copy_content(file_path, output_folder):
    with open(file_path, 'r') as file:
        content = file.read()
        endmdl_index = content.find("ENDMDL")
        endmdl_content = content[:endmdl_index + len("ENDMDL")]
        output_file_path = os.path.join(output_folder, os.path.basename(file_path))
        with open(output_file_path, 'w') as output_file:
            output_file.write(endmdl_content)
        

def extract_model1_flex(folder_name):
    input_folder = os.path.join(folder_name, "pipeline_files/2_pdbqt_out_flex")
    output_folder = os.path.join(folder_name, "pipeline_files/3_pdbqt_out_flex_m1")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file_name)
        if os.path.isfile(file_path) and file_name.endswith(".pdbqt"):
            copy_content(file_path, output_folder)

    print(f"\033[1m\033[34mExtracted Model_1 content and saved in folder: \033[91m{output_folder}\033[0m")



def process_pb_csv_flex(folder_name):
    pb_result = os.path.join(folder_name, 'pipeline_files', '2_pb_out.csv')
    pb = pd.read_csv(pb_result)
    pb = pb.drop('file', axis=1)
    pb = pb.drop(pb.index[1::2])
    pb = pb.rename(columns={'molecule': 'Name'})
    pb['passes'] = pb.iloc[:, 1:].eq('True').sum(axis=1)
    pb.to_csv(os.path.join(folder_name, 'pipeline_files', '3_pb_out.csv'), index=False)


def final_output_flex(folder_name, input_smiles, passes):
    posebusters_path = os.path.join(folder_name, 'pipeline_files', '3_pb_out.csv')
    df = pd.read_csv(posebusters_path)
    
    df1 = df[df['passes'] >= passes]

    affinity_path = os.path.join(folder_name, 'pipeline_files', '1_extract_affinity_from_pdbqt.csv')
    df2 = pd.read_csv(affinity_path)
    df2 = df2.sort_values(by='Affinity').reset_index(drop=True)

    df_temp = pd.merge(df1, df2, on='Name', how='left')

    df3 = pd.read_csv(input_smiles)

    df4 = pd.merge(df_temp[['Name']], df3[['Name', 'SMILES', 'Affinity']], on='Name', how='inner')
    df4 = df4.sort_values(by='Affinity').reset_index(drop=True)
    output_path = os.path.join(folder_name, 'output_flex.csv')
    df4.to_csv(output_path, index=False)

    filtered_count = len(df) - len(df1)
    print(f"\033[1m\033[34mCompounds filtered out by PoseBusters: \033[91m{filtered_count}\033[0m")






##############################################################################################################################
""" Process SDF file """

def process_sdf_file(sdf_file_path):
    supplier = Chem.SDMolSupplier(sdf_file_path)

    for mol in supplier:
        if mol is not None:
            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    print(f"Atom {atom.GetIdx()}: {pos.x}, {pos.y}, {pos.z}")

                img_size = (500, 500)  
                img = Draw.MolToImage(mol, size=img_size)
                img.show()
