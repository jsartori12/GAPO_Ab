
from abnumber import Chain


seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'

chain = Chain(seq, scheme='imgt', cdr_definition='imgt')

chain.cdr1_seq
chain.cdr2_seq
chain.cdr3_seq


chain.fr1_seq
chain.fr2_seq
chain.fr3_seq
chain.fr4_seq
