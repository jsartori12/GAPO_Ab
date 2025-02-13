#https://abnumber.readthedocs.io/en/latest/#position
from abnumber import Chain


#### Mudar para pegar somente a região variavel na hora de extrair os indices, para o FM nao ficar cagado. ver como selecionar somente as
#### regioes variaveis do Ab

seq = 'QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPK'

chain = Chain(seq, scheme='imgt', cdr_definition='imgt')


def extract_cdr_fm_indices(chain):
    """
    Extracts the indices of the sequences for CDRs and FMs from a Chain object.

    :param chain: An instance of the Chain class.
    :return: A dictionary with keys 'cdr_indices' and 'fm_indices', each containing a list of indices.
    """
    cdr_indices = []
    fm_indices = []

    # Iterate through all positions in the chain
    for idx, (pos, aa) in enumerate(chain.positions.items()):
        if pos.is_in_cdr():
            cdr_indices.append(idx)
        else:
            fm_indices.append(idx)

    return {
        'cdr_indices': cdr_indices,
        'fm_indices': fm_indices
    }


bla = extract_cdr_fm_indices(chain)

chain.is_light_chain()
chain.print_tall()
chain.positions()
pipi = pyrosetta.rosetta.core.pose.get_resnums_for_chain(starting_pose, "D")
''.join(starting_pose.sequence()[i-1] for i in pipi)


chain.get_region()

chain.extract_cdr_fm_indices()

for bla in chain:
    print(bla)

#### pedir pro usuario passar a cadeia pesada e leve. dai, o algoritmo extrai a sequencia ^
#### identifica as cdrs e fm, e pega os index. baseado no input, ele adiciona esses index
#### na lista de fixados para otimizados

chain.cdr1_seq
chain.cdr2_seq
chain.cdr3_seq

chain._parse_position()

chain.get_position_by_raw_index
chain.fr1_seq
chain.fr2_seq
chain.fr3_seq
chain.fr4_seq
