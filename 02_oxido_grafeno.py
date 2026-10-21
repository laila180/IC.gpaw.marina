"""
Etapa 2: Otimização do Óxido de Grafeno (GO) em Modo LCAO
Calcula o GO para as 3 supercélulas: C10, C18 e C20.
"""

from ase.io import read
from ase.optimize import BFGS
from gpaw import GPAW

def otimizar_go(arquivo_grafeno, posicoes_oxigenio, nome_saida):
    print(f"\n==========================================")
    print(f" A Iniciar Otimização GO: {nome_saida} (LCAO)")
    print(f"==========================================")
    
    # 1. Carregar a folha de grafeno otimizada
    go = read(arquivo_grafeno)
    
    # 2. Adicionar os átomos de oxigénio
    for pos in posicoes_oxigenio:
        go.append('O')
        go.positions[-1] = pos

    # 3. Configurar a calculadora GPAW em modo LCAO
    calc = GPAW(
        mode='lcao',
        basis='dzp',
        xc='PBE',
        kpts=(3, 3, 1),
        txt=f'{nome_saida}_log.txt'
    )
    go.calc = calc

    # 4. Otimizar a geometria da estrutura oxidada
    opt = BFGS(go, trajectory=f'{nome_saida}.traj')
    opt.run(fmax=0.05)

    # 5. Obter a energia e guardar a estrutura relaxada
    e_go = go.get_potential_energy()
    go.write(f'{nome_saida}_relaxed.xyz')
    print(f"-> Energia Total do {nome_saida}: {e_go:.4f} eV")
    return e_go


if __name__ == '__main__':
    # 1. GO baseado no Grafeno C10
    pos_c10 = [[1.2, 1.2, 11.4], [2.4, 2.4, 11.4]]
    e_go_c10 = otimizar_go('grafeno_C10.traj', pos_c10, 'go_C10_1O2')

    # 2. GO baseado no Grafeno C18
    pos_c18 = [[1.2, 0.7, 11.4], [2.4, 1.4, 11.4]]
    e_go_c18 = otimizar_go('grafeno_C18.traj', pos_c18, 'go_C18_1O2')

    # 3. GO baseado no Grafeno C20
    pos_c20 = [[1.2, 1.2, 11.4], [2.4, 2.4, 11.4]]
    e_go_c20 = otimizar_go('grafeno_C20.traj', pos_c20, 'go_C20_1O2')

    print("\n==========================================")
    print(" Simulações do Óxido de Grafeno concluídas com sucesso!")
    print("==========================================")
    