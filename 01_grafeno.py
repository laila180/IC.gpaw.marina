"""
Etapa 1: Otimização do Grafeno Puro (C10, C18 e C20) em Modo LCAO
Calcula a energia total e a energia por átomo para diferentes tamanhos de supercélula.
"""

from ase.build import graphene
from ase.optimize import BFGS
from gpaw import GPAW

def simular_grafeno(dimensoes, nome_sistema):
    print(f"\n==========================================")
    print(f" A Iniciar Simulação: {nome_sistema} (LCAO)")
    print(f"==========================================")
    
    # 1. Construir a supercélula de grafeno
    # A célula unitária base possui 2 átomos de carbono (C2)
    grafeno = graphene(formula='C2', a=2.46, size=dimensoes, vacuum=10.0)
    num_carbonos = len(grafeno)
    print(f"Número total de átomos de carbono no sistema: {num_carbonos}")
    
    # 2. Configurar a calculadora GPAW em modo LCAO (leve em memória RAM)
    calc = GPAW(
        mode='lcao',
        basis='dzp',
        xc='PBE',
        kpts=(3, 3, 1),
        txt=f'{nome_sistema}_log.txt'
    )
    grafeno.calc = calc
    
    # 3. Otimizar a geometria da folha
    opt = BFGS(grafeno, trajectory=f'{nome_sistema}.traj')
    opt.run(fmax=0.05)               # Critério de parada: forças < 0.05 eV/Å
    
    # 4. Obter energias e guardar a estrutura final
    e_total = grafeno.get_potential_energy()
    e_por_atomo = e_total / num_carbonos
    grafeno.write(f'{nome_sistema}_relaxed.xyz')
    
    print(f"-> Energia Total ({nome_sistema}): {e_total:.4f} eV")
    print(f"-> Energia por átomo de carbono: {e_por_atomo:.4f} eV/átomo")
    
    return e_total


if __name__ == '__main__':
    # 1. Grafeno C10 (1x5 células = 10 átomos de carbono)
    e_c10 = simular_grafeno((1, 5, 1), 'grafeno_C10')
    
    # 2. Grafeno C18 (3x3 células = 18 átomos de carbono)
    e_c18 = simular_grafeno((3, 3, 1), 'grafeno_C18')
    
    # 3. Grafeno C20 (2x5 células = 20 átomos de carbono)
    e_c20 = simular_grafeno((2, 5, 1), 'grafeno_C20')
    
    print("\n==========================================")
    print(" Simulações do grafeno puro concluídas com sucesso!")
    print("==========================================")
    