
from pathlib import Path
import pandas as pd

cur_dir = Path(__file__)
save_dir = cur_dir.parent / "result"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


from myRiboSeq import P2A2026

# P2A2026.make_fasta_snapgene(
#     snapgene_load_dir=Path('/home/ref/Reporters/plasmids'),
#     reporters=[
#         'S080-Con_one P2A No GSG',
#         'S080-Con_one P2A',
#         'S080-L_2P2A No HiBiT',
#         'S080-L_one P2A',
#         'S080-LwoHiBiT_one P2A',
#         'S080-P_one P2A',
#         'S080_pcDNA3-3xFlag-miRFP670m-P2A-cCon-P2A-EGFP-3xHA'
#     ],
#     reporter_save_names=[
#         'S080-Con_one_P2A_No_GSG',
#         'S080-Con_one_P2A',
#         'S080-L_2P2A_No_HiBiT',
#         'S080-L_one_P2A',
#         'S080-LwoHiBiT_one_P2A',
#         'S080-P_one_P2A',
#         'S080'
#     ]
# )


reporters = [
    'S080-Con_one_P2A',
    'S080','S080-LwoHiBiT_one_P2A','S080-L_2P2A_No_HiBiT', 'S080-Con_one_P2A_No_GSG',
    'S080-P_one_P2A','S080-P',
    'S080-L_one_P2A','S080-L',
    'S080-G1-L_one_P2A','S080-L-G1',
    'S080-G3-L_one_P2A','S080-L-G3',
    'S080-F_one_P2A','S080-F',
    'S080-S_one_P2A','S080-S',
    'S080-B1','S080-B2',
]

ref = P2A2026.prep_data(
    save_dir=None,
    ref_dir = f'ref/Reporters',
    data_dir="",
    reporters=reporters
)

df_smpl = pd.read_csv(save_dir / 'find_subsequence_counts' / f'exp_metadata.csv',index_col=0,header=0)

mm = 2
for rep in reporters:
    rows = df_smpl.query(f'reporter == "{rep}"')

    for fname in rows['fname'].unique():

        dataname = fname.split('_')[3]
        start_p2a = rows.query(f'fname == "{fname}"')['start_P2A'].values[0]

        for smpl in rows.query(f'fname == "{fname}"')['sample'].unique():

            colname = rows.query(f'fname == "{fname}" and sample == "{smpl}"')['column'].iloc[0]

            P2A2026.plot_indv_subsequence(
                save_dir=save_dir,
                load_query=save_dir / 'find_subsequence_query' / f'find_subsequence_{rep}_query.csv.gz',
                load_seqcnt=save_dir / 'find_subsequence_counts' / f'{fname}.csv.gz',
                ref=ref,
                reporter=rep,
                smpl=colname,
                plot_region=[start_p2a-9,start_p2a+19*3+9],
                pause_score_bg = [start_p2a, start_p2a+19*3],
                pause_score_target = [start_p2a + 18*3, start_p2a + 19*3],
                suffix=f'_P2A_{rep}_{smpl}_{dataname}',
                ylim=[-1,6]
            )