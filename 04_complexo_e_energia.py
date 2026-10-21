"""
Etapa 4: Montagem do Complexo (GO + Norfloxacino) e Cálculo da Energia de Interação
Versão com suporte a retomar cálculos existentes (Checkpoint).
"""

from os.path import exists
from ase.io import read
from ase.optimize import BFGS
from gpaw import GPAW

def simular_complexo(arquivo_go, arquivo_farma, nome_saida):
    traj_file = f'{nome_saida}.traj'
    
    print(f"\n==========================================")
    print(f" A Iniciar/Retomar Simulação: {nome_saida} (LCAO)")
    print(f"==========================================")
    
    # Verificar se já existe um ficheiro de trajetória em andamento
    if exists(traj_file):
        print(f"-> A retomar a partir do progresso anterior em {traj_file}...")
        complexo = read(traj_file, index='-1') # Carrega a última geometria calculada
    else:
        print(f"-> A iniciar o complexo de raiz...")
        go = read(arquivo_go)
        farma = read(arquivo_farma)
        farma.translate([0, 0, 3.2])
        complexo = go + farma
        complexo.center(vacuum=8.0, axis=2)
    
    # Configurar a calculadora GPAW em modo LCAO
    calc = GPAW(
        mode='lcao',
        basis='dzp',
        xc='PBE',
        kpts=(3, 3, 1),
        txt=f'{nome_saida}_log.txt'
    )
    complexo.calc = calc
    
    # Otimizar a geometria com o novo fmax (0.1 para acelerar)
    opt = BFGS(complexo, trajectory=traj_file)
    opt.run(fmax=0.1)
    
    # Obter energias finais
    e_complexo = complexo.get_potential_energy()
    go = read(arquivo_go)
    farma = read(arquivo_farma)
    e_go = go.get_potential_energy()
    e_farma = farma.get_potential_energy()
    
    e_interacao = e_complexo - (e_go + e_farma)
    complexo.write(f'{nome_saida}_relaxed.xyz')
    
    print(f"-> Energia do GO ({arquivo_go}): {e_go:.4f} eV")
    print(f"-> Energia do Fármaco:           {e_farma:.4f} eV")
    print(f"-> Energia do Complexo:          {e_complexo:.4f} eV")
    print(f"-> ENERGIA DE INTERAÇÃO (E_int): {e_interacao:.4f} eV")
    
    return e_interacao


if __name__ == '__main__':
    resultados = {}
    
    # 1. Complexo GO (C10) + Norfloxacino
    resultados['C10'] = simular_complexo('go_C10_1O2.traj', 'norfloxacino.traj', 'complexo_GO_C10_norflo')
    
    # 2. Complexo GO (C18) + Norfloxacino
    resultados['C18'] = simular_complexo('go_C18_1O2.traj', 'norfloxacino.traj', 'complexo_GO_C18_norflo')
    
    # 3. Complexo GO (C20) + Norfloxacino
    resultados['C20'] = simular_complexo('go_C20_1O2.traj', 'norfloxacino.traj', 'complexo_GO_C20_norflo')
    
    print("\n==========================================")
    print(" RESUMO FINAL - ENERGIAS DE INTERAÇÃO (LCAO)")
    print("==========================================")
    for sistema, e_int in resultados.items():
        print(f"GO_{sistema} + Norfloxacino : E_int = {e_int:.4f} eV")
    print("==========================================")
    