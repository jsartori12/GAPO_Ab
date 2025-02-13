#https://abnumber.readthedocs.io/en/latest/#position
from abnumber import Chain


seq = 'QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPK'

chain = Chain(seq, scheme='imgt', cdr_definition='imgt')


pipi = pyrosetta.rosetta.core.pose.get_resnums_for_chain(starting_pose, "D")
''.join(starting_pose.sequence()[i-1] for i in pipi)


#### pedir pro usuario passar a cadeia pesada e leve. dai, o algoritmo extrai a sequencia ^
#### identifica as cdrs e fm, e pega os index. baseado no input, ele adiciona esses index
#### na lista de fixados para otimizados

chain.cdr1_seq
chain.cdr2_seq
chain.cdr3_seq


chain.fr1_seq
chain.fr2_seq
chain.fr3_seq
chain.fr4_seq
