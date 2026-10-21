"""
Etapa 3: Otimização das Moléculas Isoladas (Gás O2 e Norfloxacino)
Versão otimizada para menor consumo de memória RAM.
"""

from ase.build import molecule
from ase.data.pubchem import pubchem_atoms_search
from ase.optimize import BFGS
from gpaw import GPAW

def otimizar_molecula(mol, nome_sistema, vacuo=3.5):
    print(f"\n====================================")
    print(f" A Iniciar Otimização: {nome_sistema}")
    print(f"====================================")

    # 1. Ajustar o vácuo ao redor da molécula
    mol.center(vacuum=vacuo)

    # 2. Configurar a calculadora GPAW em modo LCAO (baixo consumo de RAM)
    calc = GPAW(
        mode='lcao',
        basis='dzp',
        xc='PBE',
        txt=f'{nome_sistema}_log.txt'
    )
    mol.calc = calc

    # 3. Otimizar a geometria
    opt = BFGS(mol, trajectory=f'{nome_sistema}.traj')
    opt.run(fmax=0.05)

    # 4. Obter a energia total otimizada
    energia = mol.get_potential_energy()
    print(f"-> Energia Total de {nome_sistema}: {energia:.4f} eV")
    return energia


# --- Execução das Otimizações ---

# Otimização do Gás O2
print("A criar e otimizar a molécula de O2...")
gas_o2 = molecule('O2')
e_o2 = otimizar_molecula(gas_o2, 'gas_o2', vacuo=3.0)

# Otimização do Norfloxacino
print("\nA descarregar a estrutura 3D do Norfloxacino do PubChem...")
norflo = pubchem_atoms_search(name='norfloxacin')
e_norflo = otimizar_molecula(norflo, 'norfloxacino', vacuo=3.5)

print("\n====================================")
print("Otimizações concluídas com sucesso!")
print("====================================")