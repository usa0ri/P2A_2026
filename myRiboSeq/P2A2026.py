from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
import pandas as pd
import pysam
import time
import sys
import numpy as np
import re
import openpyxl
from openpyxl.styles import numbers,PatternFill
import string
from pathlib import Path, PosixPath
from tqdm import tqdm
import struct
import xmltodict
import json
from Bio import SeqIO
import gzip
import logomaker as lm
from collections import defaultdict
import codonbias as cb

import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap

sys.path.append('/home')
import myRiboSeq.myUtil as my
import myRiboSeq.myRef as myref
import myRiboSeq.myRiboBin as mybin

# import my_cython

font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
font_prop = FontProperties(fname=font_path)
plt.rcParams['font.family'] = font_prop.get_name()

dict_colors = {
    'HiBiT':'#2fbd37',
    '3xFLAG':'#ab7777',
    'start codon':'#e2e32d',
    'stop codon':'#e2e32d',
    'IGF2':'#fff42d',
    'SASP':'#92d1e3',
    'Nluc':'#a9ea81',
    'P2A':'#eacb81',
    'miRFP670':"#df8989",
    'EGFP':'#00ff00',
    '3xHA':"#74B8D8",
}

plasmid_features = [
    'lac operator',
    'M13 rev',
    'CAP binding site',
    'lac promoter',
    'AmpR promoter',
    'SV40 poly(A) signal',
    'SV40 ori',
    'CMV promoter',
    'bGH poly(A) signal',
    'SV40 promoter',
    'CMV enhancer',
    'f1 ori',
    'NeoR/KanR',
    'AmpR'
]

seq_t7 = 'TAATACGACTCACTATA'
seq_3xFLAG = 'GACTACAAGGACCACGACGGTGACTACAAGGACCACGACATCGACTACAAGGACGACGACGACAAG'
seq_HiBiT = 'GTGAGCGGCTGGCGGCTGTTCAAGAAGATTAGCTAA'

def _filter_reads(
    dict_tr:dict,
    read_len_list=[],
    softclip5 = -1,
    softclip3 = -1,
    is_full = False,
    threshold_pos_count = 1e+10
):
    n = 0;n_total = 0
    dict_tr_out = {}
    for tr,reads in tqdm(dict_tr.items()):
        n_total += len(reads)
        if is_full:
            reads = reads[
                reads['read_length'] == reads['length']
            ]
        if len(read_len_list)>0:
            reads = reads[ 
                    (reads['read_length'] >= read_len_list[0]) * (reads['read_length'] < read_len_list[1])
                ]
        
        # if (softclip3>=0) & (softclip5>=0):
        #     idx = np.array([
        #         (len(read['read_seq'])-len(read['seq'])) == (softclip5 + softclip3)
        #         for _,read in reads.iterrows()
        #     ])
        #     reads = reads.iloc[idx,:]

        if softclip5 >= 0:
            idx = np.array([
                read['read_seq'].find(read['seq']) <= softclip5
                for _,read in reads.iterrows()
            ])
            reads = reads.iloc[idx,:]
        if softclip3 >= 0:
            idx = np.array([
                (read['read_length'] - read['read_seq'].find(read['seq']) - read['length'] -1) <= softclip3
                for _,read in reads.iterrows()
            ])
            reads = reads.iloc[idx,:]
        
        reads_start = (reads['cut5'].value_counts()).to_dict()
        for k,v in reads_start.items():
            if (v > threshold_pos_count):
                reads = reads[ reads['cut5'] != k ]
                reads_start[k] = 0

        if len(reads) > 0:
            dict_tr_out[tr] = reads
            n += len(reads)

    
    print(f'{n_total-n}/{n_total} ({round((n_total-n)/n_total*100,3)}%) reads were filtered...')
    print(f'{n}/{n_total} ({round(n/n_total*100,3)}%) reads will be used...')
            
    
    return dict_tr_out

def correct_upper_lower(feature_name:str):

    keys = np.array(list(dict_colors.keys()))
    keys_upper = np.array([s.upper() for s in keys])
    if feature_name.upper() in keys_upper:
        return keys[ keys_upper == feature_name.upper() ][0]
    else:
        print(f"{feature_name} does not match")
        return ''

def save_snapgene_dict_to_dna(
        snapgene_data:dict,
        output_filepath:PosixPath
        ):
    """Save the modified SnapGene dictionary back to a .dna file."""
    title = output_filepath.stem
    with open(output_filepath, 'wb') as f:
        f.write(b'\t')  # Start byte
        f.write(struct.pack('>I', 14))  # Length
        f.write(title.encode('ascii'))  # Title
        
        # Write DNA sequence
        seq_data = snapgene_data.get("seq", "").encode('ascii')
        f.write(struct.pack('>BI', 0, len(seq_data) + 1))
        f.write(seq_data)
        
        # Write Features
        features_xml = xmltodict.unparse({"Features": {"Feature": snapgene_data["features"]}}, pretty=True)
        f.write(struct.pack('>BI', 10, len(features_xml)))
        f.write(features_xml.encode('utf-8'))

def add_feature_to_snapgene(
        snapgene_data:dict, 
        name:str, 
        start:int, 
        end:int, 
        strand='+', 
        feature_type='misc_feature', 
        color='#FF0000'
    ):
    """Add a new feature to the SnapGene data dictionary."""
    new_feature = {
        "start": start, 
        "end": end,
        "strand": strand,
        "type": feature_type,
        "name": name,
        "color": color,
        "textColor": "black",
        "segments": [{"@range": f"{start}-{end}", "@color": color}],
        "row": 0,
        "isOrf": False,
        "qualifiers": {"label": name, "note": ["color: " + color]}
    }
    snapgene_data["features"].append(new_feature)
    return snapgene_data

class Ref:

    def __init__(
        self,
        data_dir:PosixPath,
        ref_dir:PosixPath,
        reporters:list,
        aligned_seq_path:PosixPath = Path()
        ):

        data_dir = Path(data_dir)
        ref_dir = Path(ref_dir)

        if data_dir != Path():
            self.data_dir = data_dir
            exp_metadata_file = data_dir / "exp_metadata.csv"
            '''parse metadata of the experiment'''
            self.exp_metadata = myref.ExpMetadata(exp_metadata_file)
            print("experimental metadata has been loaded...")
        
        self.ref_reporters = {}
        self.seq_reporters = {}
        for r in reporters:
            with open(ref_dir / 'snapgene_json' / (r + '.json'),'r') as f:
                self.ref_reporters[r] = json.load(f)
            with open(ref_dir / 'fasta' / (r + '.fa'),'r') as f:
                for rec in SeqIO.parse(f,'fasta'):
                    self.seq_reporters[r] = str(rec.seq)
        
        self.dict_id = dict(zip(
            reporters,
            list(range(len(reporters)))
        ))

        if aligned_seq_path != Path():
            dict_align = {}
            with open(aligned_seq_path,'r') as f:
                for rec in SeqIO.parse(f,'fasta'):
                    dict_align[rec.id] = str(rec.seq)
        
            self.dict_align = dict_align

  

def prep_data(
    save_dir:Path,
    data_dir:Path,
    ref_dir:Path,
    reporters:list,
    aligned_seq_path=Path()
):
    ref = Ref(
        data_dir=data_dir,
        ref_dir=ref_dir,
        reporters=reporters,
        aligned_seq_path=aligned_seq_path
    )

    if save_dir is None:
        return ref

    save_dir = save_dir / 'prep_data'
    if not save_dir.exists():
        save_dir.mkdir()

    n_reads_list = []
    for a in ref.exp_metadata.df_metadata['align_files']:
        smpl_name = ref.exp_metadata.df_metadata.query(f'align_files == "{a}"').sample_name.iloc[-1]
        print(f'preparing data for {smpl_name}...')
        infile = pysam.AlignmentFile(ref.data_dir / a)

        obj = mybin.myBinRiboReporter(
            data_dir=data_dir,
            smpl=smpl_name,
            save_dir=save_dir,
            dict_id=ref.dict_id
        )
        obj.encode(infile,ref)
        n_reads_list.append(obj.n_reads)

    pd.DataFrame(n_reads_list,index=ref.exp_metadata.df_metadata['sample_name'],columns=['total_reads']).\
        to_csv(save_dir / 'summary.csv.gz')
    
    return ref

def indiv_plot(
    save_dir,
    load_dir,
    fname,
    smpls,
    ref,
    thres_pos=0,
    ylim_now=[],
    is_full=False,
    softclip5=-1,
    softclip3=-1,
    offset=0
):
    save_dir = save_dir / f'indiv_plot'
    if not save_dir.exists():
        save_dir.mkdir()

    num_smpl = len(smpls)

    outfile_name = save_dir / f'indiv_plot_{fname}'
    pdf = PdfPages(outfile_name.with_suffix(".pdf"))   

    dict_tr_list = {};rep_names = []
    for i,s in enumerate(smpls):
        obj = mybin.myBinRiboReporter(
            data_dir=ref.data_dir,
            smpl=s,
            save_dir=load_dir,
            dict_id=ref.dict_id
        )
        obj.decode()
        df_data,dict_tr = obj.make_df(is_seq=True)
        dict_tr = _filter_reads(
            dict_tr=dict_tr,
            softclip5=softclip5,
            softclip3=softclip3,
            is_full=is_full,
            threshold_pos_count=1e+10
        )

        df_data = []
        for k,v in dict_tr.items():
            if len(df_data) == 0:
                df_data = v
            else:
                df_data = pd.concat(df_data,v,axis=0)
            start_mrna = 0
            if "3'UTR" in ref.ref_reporters[k].keys():
                end_mrna = ref.ref_reporters[k]["3'UTR"]['end']
            elif "3' UTR" in ref.ref_reporters[k].keys():
                end_mrna = ref.ref_reporters[k]["3' UTR"]['end']
            else:
                end_mrna = len(ref.ref_reporters[k]['seq'])
            df_data = df_data[ (df_data['cut5']>start_mrna-offset)*(df_data['cut5']<end_mrna-offset) ]
            df_data['cut5'] -= start_mrna
            df_data['cut3'] -= start_mrna
        print(f'{len(df_data)} reads of {s} are used for plotting...')
        dict_tr_list[s] = df_data
        rep_names += list(np.unique(list(dict_tr.keys())))
        if len(df_data)>0:
            df_data.to_csv(save_dir / f'df_data_{s}_{fname}.csv.gz')
    
    rep_name = np.unique(rep_names)[0]

    if num_smpl>1:

        max_len = np.amax([
            len(ref.ref_reporters[df_data['tr_id'].iloc[0]]['seq'])
            for df_data in dict_tr_list.values()
        ])
        if max_len > 700:
            w = 18
        else:
            w = 12
        fig, ax = plt.subplots(num_smpl,1,figsize=(w,1.5*num_smpl),sharex=True,sharey=True)
    
        ax[-1].set_xlabel(f'Position along the {rep_name} (nt)') 
        ax[-1].set_ylabel('Read counts')   

        for i,(s,df_data) in enumerate(dict_tr_list.items()):

            ax[i].set_title(s)

            if len(df_data)==0:
                continue

            df_cnt = df_data.groupby(['cut5']).apply(len)
            df_cnt = df_cnt.set_axis( df_cnt.index + offset, axis=0)
            if thres_pos > 0:
                df_cnt[df_cnt>thres_pos] = df_cnt[ df_cnt<=thres_pos ].max(axis=0)
            tr = df_data['tr_id'].unique()[0]
            
            for region_dict in ref.ref_reporters[tr]['features']:
                region_name = region_dict['name']
                if 'rich' in region_name:
                    color = '#2deade'
                ax[i].axvspan(region_dict['start'],region_dict['end'],color=color)
        
            idx_frame0 = (df_cnt.index - (ref.ref_reporters[tr]['start codon']['start'])) % 3 == 0
            ax[i].vlines(df_cnt.index[~idx_frame0],0,df_cnt.values[~idx_frame0],colors='#808080',lw=1)
            ax[i].vlines(df_cnt.index[idx_frame0],0,df_cnt.values[idx_frame0],colors='#FF0000',lw=1)
            
            if len(ylim_now)>0:
                ax[i].set_ylim(ylim_now[0],ylim_now[1])
            
            aa_seq = my._seq2aa((ref.ref_reporters[tr]['seq'][
                ref.ref_reporters[tr]['start codon']['start'] : ref.ref_reporters[tr]['stop codon']['end']
            ]).upper())
            for j,aa in enumerate(aa_seq):
                ax[i].text(
                    j*3+1 + ref.ref_reporters[tr]['start codon']['start'], 0, aa,
                    va='top', ha='center',fontsize=5
                )
            
            df_cnt.to_csv(save_dir / f'indiv_plot_{fname}_{s}.csv.gz')
        
        j = 0
        for region_dict in ref.ref_reporters[tr]['features']:
            ax[i].text(
                (region_dict['start'] + region_dict['end']) / 2,
                ax[i].get_ylim()[1]*0.75 + (-1)**(j%2)*ax[i].get_ylim()[1]*0.1,
                region_name,
                ha='center'
            )
            j += 1

        fig.tight_layout()
        # plt.show()
        fig.savefig(pdf, format='pdf')
        fig.savefig(outfile_name.with_suffix('.png'),dpi=300)
        plt.close()
        pdf.close()

        
    else:

        tr = df_data['tr_id'].unique()[0]
        
        if len(ref.seq_reporters[tr]) > 700:
            w = 18
        else:
            w = 12
        fig, ax = plt.subplots(1,1,figsize=(w,1.8),sharex=True,sharey=True)

        ax.set_xlabel(f'Position along the {rep_name} (nt)') 
        ax.set_ylabel('Read counts')

        df_cnt = df_data.groupby(['cut5']).apply(len)
        df_cnt = df_cnt.set_axis( df_cnt.index + offset, axis=0)
        if thres_pos > 0:
            df_cnt[df_cnt>thres_pos] = df_cnt[ df_cnt<=thres_pos ].max(axis=0)
            
        for region_dict in ref.ref_reporters[tr]['features']:
            region_name = region_dict['name']
            if 'rich' in region_name:
                color = '#2deade'
            ax.axvspan(region_dict['start'],region_dict['end'],color=region_dict['color'])

            if region_name == 'start codon':
                pos_start_codon = region_dict['start']
            if region_name == 'stop codon':
                pos_stop_codon = region_dict['start']
        
        idx_frame0 = ((df_cnt.index - pos_start_codon) % 3 == 0) *\
            (df_cnt.index >= pos_start_codon) *\
            (df_cnt.index <= pos_stop_codon)
        ax.vlines(df_cnt.index[~idx_frame0],0,df_cnt.values[~idx_frame0],colors='#808080',lw=1)
        ax.vlines(df_cnt.index[idx_frame0],0,df_cnt.values[idx_frame0],colors='#FF0000',lw=1)
        ax.set_title(s)

        if len(ylim_now)>0:
            ax.set_ylim(ylim_now[0],ylim_now[1])
            
        aa_seq = my._seq2aa((ref.seq_reporters[tr][
            pos_start_codon: pos_stop_codon
        ]).upper())
        fontsize_now = 5 if len(aa_seq)<100 else 3
        for j,aa in enumerate(aa_seq):
            ax.text(
                j*3+1 + pos_start_codon, 0, aa,
                va='top', ha='center',fontsize=fontsize_now
            )
            
        j = 0
        for region_dict in ref.ref_reporters[tr]['features']:
            ax.text(
                (region_dict['start'] + region_dict['end']) / 2,
                ax.get_ylim()[1]*0.75 + (-1)**(j%2)*ax.get_ylim()[1]*0.1,
                region_dict['name'],
                ha='center'
            )
            j += 1
        
        fig.tight_layout()
        # plt.show()
        fig.savefig(pdf, format='pdf')
        fig.savefig(outfile_name.with_suffix('.png'),dpi=300)
        plt.close()
        pdf.close()

        df_cnt.to_csv(outfile_name.with_suffix('.csv.gz'))


def indiv_plot_edit(
    save_dir,
    dict_load_files,
    outfile_name,
    tr:str,
    ref,
    ylim_now=[]
):
    save_dir = save_dir / f'indiv_plot'
    if not save_dir.exists():
        save_dir.mkdir()
    
    num_smpl = len(dict_load_files)

    start_mrna = 0
    end_mrna = len(ref.ref_reporters[tr]['seq'])
    
    pdf = PdfPages(outfile_name)
    fig, ax = plt.subplots(num_smpl,1,figsize=(10,1.5*num_smpl),sharex=True,sharey=True)
    ax[-1].set_xlabel(f'Position along the {tr} (nt)') 
    ax[-1].set_ylabel('Read counts')   


    for i,(s,load_file) in enumerate(dict_load_files.items()):
        df_cnt = pd.read_csv(load_file,index_col=0,header=0)

        ax[i].set_title(s)

        for region_name,region_dict in ref.ref_reporters[tr].items():
            
            if region_name == 'seq':
                continue
            if (region_dict['start'] > start_mrna) * (region_dict['end'] < end_mrna):
                color = dict_colors.get(region_name,region_dict['color'])
                if 'rich' in region_name:
                    color = '#2deade'
                ax[i].axvspan(region_dict['start']-start_mrna,region_dict['end']-start_mrna,color=color)
    
        idx_frame0 = (df_cnt.index - (ref.ref_reporters[tr]['start codon']['start']-start_mrna)) % 3 == 0
        ax[i].vlines(df_cnt.index[~idx_frame0],0,df_cnt.values[~idx_frame0],colors='#808080',lw=1)
        ax[i].vlines(df_cnt.index[idx_frame0],0,df_cnt.values[idx_frame0],colors='#FF0000',lw=1)
        
        if len(ylim_now)>0:
            ax[i].set_ylim(ylim_now[0],ylim_now[1])
        
        aa_seq = my._seq2aa((ref.ref_reporters[tr]['seq'][
            ref.ref_reporters[tr]['start codon']['start'] : ref.ref_reporters[tr]['stop codon']['end']
        ]).upper())
        for j,aa in enumerate(aa_seq):
            ax[i].text(
                j*3+1 + ref.ref_reporters[tr]['start codon']['start'], 0, aa,
                va='top', ha='center',fontsize=5
            )
    
    j = 0
    for region_name,region_dict in ref.ref_reporters[tr].items():
        if region_name == 'seq':
            continue
        if (region_dict['start'] > start_mrna) * (region_dict['end'] < end_mrna):
            ax[i].text(
                (region_dict['start']-start_mrna + region_dict['end']-start_mrna) / 2,
                ax[i].get_ylim()[1]*0.75 + (-1)**(j%2)*ax[i].get_ylim()[1]*0.1,
                region_name,
                ha='center'
            )
            j += 1

    fig.tight_layout()
    plt.show()
    fig.savefig(pdf, format='pdf')
    fig.savefig(outfile_name.with_suffix('.png'),dpi=300)
    plt.close()
    pdf.close()

def make_fasta_AA_snapgene(
    snapgene_load_dir:PosixPath,
    reporters:list,
    reporter_save_names:list,
    fasta_dir=Path('/home/ref/Reporters/fasta'),
    snapgene_dir=Path('/home/ref/Reporters/snapgene'),
    snapgene_json_dir=Path('/home/ref/Reporters/snapgene_json')
):

    assert len(reporters) == len(reporter_save_names),\
        "reporters and reporter_save_names should have the same length"

    for reporter, reporter_save_name in zip(reporters,reporter_save_names):
        dict_snapgene = snapgene_file_to_dict(snapgene_load_dir / (reporter + '.dna'))

        # see whether mRNA region is defined or not
        n_mrna_feature = np.sum([ ft['name'] in ['mRNA',reporter_save_name] for ft in dict_snapgene['features'] ])
        if n_mrna_feature > 1:
            # error
            raise Exception(
                f"mRNA feature is already defined for {reporter_save_name}. "
                "Please remove it before running this function."
            )
        is_mrna_feature = n_mrna_feature > 0
        if is_mrna_feature:
            for dict_feature in dict_snapgene['features']:
                if dict_feature['name'] in ['mRNA',reporter_save_name]:
                    seq_new = dict_snapgene['seq'][ dict_feature['start'] : dict_feature['end'] ]
                    start_mrna = dict_feature['start']
                    end_mrna = dict_feature['end']
                    is_mrna_feature = True
        else:
            # detect XbaI site (mRNA should be end just before XbaI)
            pos_xbai = dict_snapgene['seq'].upper().find('TCTAGA')
            if pos_xbai == -1:
                print("no XbaI site")
                end_mrna = dict_snapgene['dna']['length']
            start_mrna = 0
            end_mrna = pos_xbai-1
            seq_new = dict_snapgene['seq'][start_mrna:end_mrna]
            
        # trim T7 promoter sequence
        seq_new = seq_new.upper()
        is_t7 = seq_t7 in seq_new
        if is_t7:
            seq_new = seq_new.split(seq_t7)[-1]
        start_mrna = dict_snapgene['seq'].upper().find(seq_new)
        
        # detect kozak + start codon
        pos_start_codon = seq_new.find('ACCATG') + 3
        if pos_start_codon == 2:
            print("no start codon")
        
        features_list = [];pos_max = 0
        if seq_new != '':
            # detect start/stop codons
            # make new feature list
            for dict_feature in dict_snapgene['features']:
                if dict_feature['name'] in plasmid_features:
                    continue
                if dict_feature['name'] in ['start codon','stop codon']:
                    continue
                if (dict_feature['start'] >= pos_start_codon+start_mrna) * (dict_feature['end'] <= end_mrna):
                    name_correct = dict_feature['name']
                    if name_correct not in dict_colors.keys():
                        name_correct = correct_upper_lower(feature_name=dict_feature['name'])
                    if name_correct == '':
                        continue
                    dict_feature['color'] = dict_colors[name_correct]
                    dict_feature['name'] = name_correct
                    if dict_feature['name'] == '3xFLAG':
                        start_pos = seq_new.find(seq_3xFLAG)
                        end_pos = start_pos + len(seq_3xFLAG)
                    elif dict_feature['name'] == 'HiBiT':
                        start_pos = seq_new.find(seq_HiBiT)
                        end_pos = start_pos + len(seq_HiBiT)
                    else:
                        start_pos = dict_feature['start'] - start_mrna
                        end_pos = dict_feature['end'] - start_mrna
                    dict_feature['start'] = start_pos
                    dict_feature['end'] = end_pos
                    features_list.append(dict_feature)
                    pos_max = int(np.amax([pos_max,dict_feature['end']]))
        else:
            print("no mRNA features")
        dict_snapgene_new = dict_snapgene.copy()
        dict_snapgene_new['features'] = features_list
        dict_snapgene_new['seq'] = seq_new
        dict_snapgene_new['dna']['length'] = len(seq_new)

        assert 'stop codon' not in dict_snapgene_new['features']
        assert 'start codon' not in dict_snapgene_new['features']

        # add start/stop codon features
        dict_snapgene_new = add_feature_to_snapgene(
            snapgene_data=dict_snapgene_new,
            name='start codon',
            start=pos_start_codon,
            end=pos_start_codon+3,
            color=dict_colors['start codon'])
        
        aa_tmp = my._seq2aa(seq_new[pos_start_codon:])
        pos_stop_codon = int(np.where(aa_tmp == '_')[0][0]*3 + pos_start_codon)
        assert seq_new[pos_start_codon:pos_start_codon+3] == 'ATG'
        assert seq_new[pos_stop_codon:pos_stop_codon+3] == 'TAA'

        dict_snapgene_new = add_feature_to_snapgene(
            snapgene_data=dict_snapgene_new,
            name='stop codon',
            start=pos_stop_codon,
            end=pos_stop_codon +3,
            color=dict_colors['stop codon'])

        with open(snapgene_json_dir / (reporter_save_name + '.json'), 'w') as f:
            json.dump(dict_snapgene_new,f)
        # save_snapgene_dict_to_dna(dict_snapgene_new,snapgene_dir / (reporter_save_name + '.dna'))

        outfa_name = fasta_dir / (reporter_save_name + '.fa')
        with open(outfa_name,'w') as f:
            f.write(f'>{reporter_save_name}\n')
            f.write(dict_snapgene_new['seq'])
        
        print(f'fasta and json files were saved for {reporter_save_name}')


def make_fasta_snapgene(
    snapgene_load_dir:PosixPath,
    reporters:list,
    reporter_save_names:list,
    fasta_dir=Path('/home/ref/Reporters/fasta'),
    snapgene_dir=Path('/home/ref/Reporters/snapgene'),
    snapgene_json_dir=Path('/home/ref/Reporters/snapgene_json')
):

    assert len(reporters) == len(reporter_save_names),\
        "reporters and reporter_save_names should have the same length"

    for reporter, reporter_save_name in zip(reporters,reporter_save_names):
        dict_snapgene = snapgene_file_to_dict(snapgene_load_dir / (reporter + '.dna'))

        # trim T7 promoter sequence
        dict_snapgene['seq'] = dict_snapgene['seq'].upper()
        seq_new = dict_snapgene['seq']
        is_t7 = seq_t7 in seq_new
        if is_t7:
            seq_new = seq_new.split(seq_t7)[-1]
        start_mrna = dict_snapgene['seq'].upper().find(seq_new)

        # end_mrna is before bGH poly(A) signal
        if 'bGH poly(A) signal' not in [ ft['name'] for ft in dict_snapgene['features'] ]:
            print("no bGH poly(A) signal")
            end_mrna = dict_snapgene['dna']['length']
        else:
            end_mrna = [ ft['start'] for ft in dict_snapgene['features'] if ft['name'] == 'bGH poly(A) signal' ][0]-1
        # try to detect features called 3'UTR' or 'mRNA'
        for dict_feature in dict_snapgene['features']:
            if '3\'UTR' in dict_feature['name'] or 'mRNA' in dict_feature['name']:
                end_mrna = dict_feature['end'] -1
        seq_new = dict_snapgene['seq'][start_mrna:end_mrna]

        # detect kozak + start codon
        pos_start_codon = seq_new.find('ACCATG') + 3
        if pos_start_codon == 2:
            print("no start codon")
            for ft in dict_snapgene['features']:
                if ft['name'] == 'start codon':
                    pos_start_codon = ft['start'] - start_mrna
        
        features_list = [];pos_max = 0
        if seq_new != '':
            # detect start/stop codons
            # make new feature list
            for dict_feature in dict_snapgene['features']:
                if dict_feature['name'] in plasmid_features:
                    continue
                if dict_feature['name'] in ['start codon','stop codon']:
                    continue
                if (dict_feature['start'] >= pos_start_codon+start_mrna-3) * (dict_feature['end'] <= end_mrna):
                    name_correct = dict_feature['name']
                    if name_correct not in dict_colors.keys():
                        name_correct = correct_upper_lower(feature_name=dict_feature['name'])
                    if name_correct == '':
                        continue
                    dict_feature['color'] = dict_colors[name_correct]
                    dict_feature['name'] = name_correct
                    start_pos = dict_feature['start'] - start_mrna
                    end_pos = dict_feature['end'] - start_mrna
                    dict_feature['start'] = start_pos
                    dict_feature['end'] = end_pos
                    features_list.append(dict_feature)
                    pos_max = int(np.amax([pos_max,dict_feature['end']]))
        else:
            print("no mRNA features")
        dict_snapgene_new = dict_snapgene.copy()
        dict_snapgene_new['features'] = features_list
        dict_snapgene_new['seq'] = seq_new
        dict_snapgene_new['dna']['length'] = len(seq_new)

        assert 'stop codon' not in dict_snapgene_new['features']
        assert 'start codon' not in dict_snapgene_new['features']

        # add start/stop codon features
        dict_snapgene_new = add_feature_to_snapgene(
            snapgene_data=dict_snapgene_new,
            name='start codon',
            start=pos_start_codon,
            end=pos_start_codon+3,
            color=dict_colors['start codon'])
        
        aa_tmp = my._seq2aa(seq_new[pos_start_codon:])
        pos_stop_codon = int(np.where(aa_tmp == '_')[0][0]*3 + pos_start_codon)
        assert seq_new[pos_start_codon:pos_start_codon+3] == 'ATG'
        assert seq_new[pos_stop_codon:pos_stop_codon+3] in ['TAG','TAA','TGA']

        dict_snapgene_new = add_feature_to_snapgene(
            snapgene_data=dict_snapgene_new,
            name='stop codon',
            start=pos_stop_codon,
            end=pos_stop_codon +3,
            color=dict_colors['stop codon'])

        with open(snapgene_json_dir / (reporter_save_name + '.json'), 'w') as f:
            json.dump(dict_snapgene_new,f)
        # save_snapgene_dict_to_dna(dict_snapgene_new,snapgene_dir / (reporter_save_name + '.dna'))

        outfa_name = fasta_dir / (reporter_save_name + '.fa')
        with open(outfa_name,'w') as f:
            f.write(f'>{reporter_save_name}\n')
            f.write(dict_snapgene_new['seq'])
        
        print(f'fasta and json files were saved for {reporter_save_name}')


def make_suffix_array(
    save_dir:PosixPath,
    ref:Ref
):
    from pydivsufsort import divsufsort, kasai

    for rep_name,seq in ref.seq_reporters.items():
        sa = divsufsort(seq)
        lcp = kasai(seq,sa)
        np.savez_compressed(save_dir / (rep_name+'.npz'),sa=sa,lcp=lcp)

def find_subsequence(
    save_dir:PosixPath,
    dict_fastq:dict,
    reporter_name:str,
    suffix='',
    k=29
):
    save_dir = save_dir / 'find_subsequence'
    if not save_dir.exists():
        save_dir.mkdir()
    
    with open(f'/home/ref/Reporters/snapgene_json/{reporter_name}.json','r') as f:
        dict_tmp = json.load(f)
    
    dict_seq_query = {}
    for i in range(len(dict_tmp['seq'])-k):
        dict_seq_query[ f'{reporter_name}_{i}' ] =dict_tmp['seq'][i : i+k]
    pd.DataFrame().from_dict(dict_seq_query,orient='index').to_csv(save_dir / f'find_subsequence_{suffix}_query.csv.gz')
    
    df = pd.DataFrame()

    for smpl,load_path in dict_fastq.items():

        dict_cnt = dict.fromkeys(dict_seq_query.keys(),0)

        seqs = []
        with gzip.open(load_path,'rt') as f:
            for j,rec in tqdm(enumerate(SeqIO.parse(f,'fastq')),desc='loading fastq data'):
                seqs.append(str(rec.seq))
        
        seqs_uniq,cnt_uniq = np.unique(seqs,return_counts=True)

        seq_list = {}
        for seq_name in dict_seq_query.keys():
            seq_list[seq_name] = []

        for seq,cnt in tqdm(zip(seqs_uniq,cnt_uniq),desc='searching query sequences from fastq data'):
            seq_names_list = []
            for seq_name,seq_query in dict_seq_query.items():
                if seq_query in seq:
                    seq_names_list.append(seq_name)
            if len(seq_names_list) == 0:
                continue
            elif len(seq_names_list) == 1:
                seq_name = seq_names_list[0]            
            else:
                pos_list = np.array([ seq.find(dict_seq_query[s]) for s in seq_names_list ])
                idx = pos_list==3
                if np.any(idx):
                    seq_name = np.array(seq_names_list)[ idx ][0]
                else:
                    seq_name = seq_names_list[0]
            dict_cnt[seq_name] += cnt
            seq_list[seq_name].append(seq)
        
        df = pd.merge(
            df,pd.DataFrame().from_dict(dict_cnt,orient='index').set_axis([smpl],axis=1),
            how='outer',left_index=True,right_index=True
        )

        with gzip.open(save_dir / f'find_subsequence_{smpl}_{suffix}.json.gz','wt') as f:
            json.dump(seq_list,f)

       
    outfile_name = save_dir / f'find_subsequence_{suffix}.csv.gz'
    df.to_csv(outfile_name)
    pd.DataFrame().from_dict(dict_seq_query,orient='index').to_csv(save_dir / f'find_subsequence_{suffix}_query.csv.gz')


def find_subsequence_edlib(
    save_dir:PosixPath,
    dict_fastq:dict,
    reporter_name:str,
    suffix='',
    k=29,
    mm=1,
    is_allow_multimap = False
):
    save_dir = save_dir / 'find_subsequence_edlib'
    if not save_dir.exists():
        save_dir.mkdir()
    
    with open(f'/home/ref/Reporters/snapgene_json/{reporter_name}.json','r') as f:
        dict_tmp = json.load(f)
    seq_ref = dict_tmp['seq']

    dict_seq_query = {
        f'{reporter_name}_{i}':seq_ref[i:i+k]
        for i in range(len(seq_ref)-k+1)
    }
    pd.DataFrame().from_dict(dict_seq_query,orient='index').to_csv(save_dir / f'find_subsequence_{suffix}_query.csv.gz')

    df = pd.DataFrame()

    for smpl,load_path in dict_fastq.items():

        print('generating k-mers from reads...')
        kmer_dict, read_dict = my.run_cython_fastq(
            load_path=load_path,
            myfun=lambda pbar: my_cython.extract_kmers(load_path.as_posix().encode('utf-8'),k,pbar)
        )

        # kmer_dict = defaultdict(list)
        # with gzip.open(load_path,'rt') as f:
        #     for j,rec in tqdm(enumerate(SeqIO.parse(f,'fastq')),desc='loading fastq data & making k-mers'):
        #         read = str(rec.seq)
        #         l = len(read)
        #         for i in range(l-k+1):
        #             kmer = read[i:i+k]
        #             kmer_dict[kmer].append((rec.id, i))

        print('searching k-mers from reference...')
        kmer_matches, readid_matches = my.run_cython_list(
            input_list=kmer_dict,
            myfun = lambda pbar: my_cython.search_kmers(kmer_dict,seq_ref,mm,pbar)
        )

        # kmer_matches = defaultdict(list)
        # for kmer in tqdm(kmer_dict.keys(),desc='searching query sequences from fastq data'):
        #     res = edlib.align(kmer,seq_ref,mode="HW",task="location")
        #     if res['editDistance'] > mm:
        #         continue
        #     for (start,end) in res['locations']:
        #         kmer_matches[kmer].append([end+1-k,end+1,res['editDistance']])
                # select only 1 alignment per read
        dict_cnt = {mm_:0 for mm_ in range(mm+1)}
        kmer_matches_out = defaultdict(list)

        if is_allow_multimap:
            for readid,alignments in readid_matches.items():
                for v in alignments:
                    kmer_matches_out[f'{reporter_name}_{v[1]}'].append(read_dict[readid])
        else:
            for readid,alignments in readid_matches.items():
                if len(alignments)==1:
                    continue
                # first, select edit_distance == 0
                idx_keep = np.array([v[-1] for v in alignments]) == 0
                if np.sum(idx_keep) == 0:
                    # if no alignment is edit_distance == 0, select kmer_start_pos == 3
                    idx_keep = np.array([v[0] for v in alignments]) == 3
                    # if no alignment is kmer_start_pos == 3, select the first alignment
                    if np.sum(idx_keep) == 1:
                        idx_keep = idx_keep
                    else:
                        idx_keep = np.array([True]+[False]*(len(alignments)-1))
                    
                elif np.sum(idx_keep) == 1:
                    idx_keep = idx_keep 
                else:
                    alignments = [v for v,i in zip(alignments,idx_keep) if i]
                    idx_keep = np.array([v[0] for v in alignments]) == 3
                    # if no alignment is kmer_start_pos == 3, select the first alignment
                    if np.sum(idx_keep) == 1:
                        idx_keep = idx_keep
                    else:
                        idx_keep = np.array([True]+[False]*(len(alignments)-1))

                for v,is_keep in zip(alignments,idx_keep):
                    if is_keep:
                        kmer_matches_out[f'{reporter_name}_{v[1]}'].append(read_dict[readid])
                        dict_cnt[v[-1]] += 1  
                        continue
                    if v[1] < 0:
                        continue
            
        pos_cnt = {i:0 for i in range(len(seq_ref)-k+1)}
        for pos_str,seqs in kmer_matches_out.items():
            pos_cnt[ int(pos_str.split('_')[-1]) ] = len(seqs)

        df = pd.merge(
            df,pd.DataFrame().from_dict(pos_cnt,orient='index').set_axis([smpl],axis=1),
            how='outer',left_index=True,right_index=True
        )
        pd.DataFrame().from_dict(dict_cnt,orient='index').to_csv(save_dir / f'find_subsequence_mm{mm}_{smpl}_{suffix}_cnt.csv.gz')
        
        outfile_name = save_dir / f'find_subsequence_mm{mm}_{smpl}_{suffix}.json.gz'
        with gzip.open(outfile_name,'wt') as f:
            json.dump(kmer_matches_out,f)
    
    outfile_name = save_dir / f'find_subsequence_mm{mm}_{suffix}.csv.gz'
    df.to_csv(outfile_name)


def find_subsequence_edlib_iter(
    save_dir:PosixPath,
    dict_fastq:dict,
    reporter_name:str,
    suffix='',
    k=29,
    mm=1
):
    import edlib

    save_dir = save_dir / 'find_subsequence_edlib'
    if not save_dir.exists():
        save_dir.mkdir()
    
    if Path(reporter_name).exists():
        if reporter_name.endswith('.fa'):
            with open(reporter_name,'r') as f:
                lines = f.readlines()
                seq_ref = ''.join([ line.strip() for line in lines if not line.startswith('>') ])
    else:
        with open(f'/home/ref/Reporters/snapgene_json/{reporter_name}.json','r') as f:
            dict_tmp = json.load(f)
        seq_ref = dict_tmp['seq']

    dict_seq_query = {
        f'{reporter_name}_{i}':seq_ref[i:i+k]
        for i in range(len(seq_ref)-k+1)
    }
    pd.DataFrame().from_dict(dict_seq_query,orient='index').to_csv(save_dir / f'find_subsequence_{suffix}_query.csv.gz')

    df = pd.DataFrame()

    for smpl,load_path in dict_fastq.items():

        kmer_matches = defaultdict(list);readid_matches = defaultdict(list)
        read_dict = {}
        with gzip.open(load_path,'rt') as f:
            for j,rec in tqdm(enumerate(SeqIO.parse(f,'fastq')),desc='loading fastq data & making k-mers & finding k-mers from reference'):
                read = str(rec.seq)
                l = len(read)
                for i in range(l-k+1):
                    kmer = read[i:i+k]
                    res = edlib.align(kmer,seq_ref,mode="HW",task="location")
                    if res['editDistance'] > mm:
                        continue
                    read_dict[rec.id] = read
                    for (start,end) in res['locations']:
                        kmer_matches[kmer].append([end+1-k,end+1,res['editDistance']])
                        readid_matches[rec.id].append([i,end+1-k,end+1,res['editDistance']])

        # select only 1 alignment per read
        dict_cnt = {mm_:0 for mm_ in range(mm+1)}
        kmer_matches_out = defaultdict(list)
        for readid,alignments in readid_matches.items():
            if len(alignments)==1:
                continue
            # first, select edit_distance == 0
            idx_keep = np.array([v[-1] for v in alignments]) == 0
            if np.sum(idx_keep) == 0:
                # if no alignment is edit_distance == 0, select kmer_start_pos == 3
                idx_keep = np.array([v[0] for v in alignments]) == 3
                # if no alignment is kmer_start_pos == 3, select the first alignment
                if np.sum(idx_keep) == 1:
                    idx_keep = idx_keep
                else:
                    idx_keep = np.array([True]+[False]*(len(alignments)-1))
                
            elif np.sum(idx_keep) == 1:
                idx_keep = idx_keep 
            else:
                alignments = [v for v,i in zip(alignments,idx_keep) if i]
                idx_keep = np.array([v[0] for v in alignments]) == 3
                # if no alignment is kmer_start_pos == 3, select the first alignment
                if np.sum(idx_keep) == 1:
                    idx_keep = idx_keep
                else:
                    idx_keep = np.array([True]+[False]*(len(alignments)-1))

            for v,is_keep in zip(alignments,idx_keep):
                if is_keep:
                    kmer_matches_out[f'{reporter_name}_{v[1]}'].append(read_dict[readid])
                    dict_cnt[v[-1]] += 1  
                    continue
                if v[1] < 0:
                    continue
        
        pos_cnt = {i:0 for i in range(len(seq_ref)-k+1)}
        for pos_str,seqs in kmer_matches_out.items():
            pos_cnt[ int(pos_str.split('_')[-1]) ] = len(seqs)

        df = pd.merge(
            df,pd.DataFrame().from_dict(pos_cnt,orient='index').set_axis([smpl],axis=1),
            how='outer',left_index=True,right_index=True
        )
        pd.DataFrame().from_dict(dict_cnt,orient='index').to_csv(save_dir / f'find_subsequence_mm{mm}_{smpl}_{suffix}_cnt.csv.gz')
        
        outfile_name = save_dir / f'find_subsequence_mm{mm}_{smpl}_{suffix}.json.gz'
        with gzip.open(outfile_name,'wt') as f:
            json.dump(kmer_matches_out,f)
    
    outfile_name = save_dir / f'find_subsequence_mm{mm}_{suffix}.csv.gz'
    df.to_csv(outfile_name)

def find_poly_anywhere(seq: str, base='A', min_len=10, min_frac=0.85):
    """
    return the first region in the sequence that has length >= min_len and base fraction >= min_frac.
    if not found, return None. O(n^2) but should be fast enough for ~150bp sequences.
    """
    s = seq.upper(); t = base.upper(); n = len(s)
    for L in range(0, n - min_len + 1):
        a = 0; best_R = None
        for R in range(L, n):
            if s[R] == t: a += 1
            win_len = R - L + 1
            if win_len >= min_len and (a / win_len) >= min_frac:
                best_R = R + 1  # 右端は非包含
        if best_R is not None:
            return (L, best_R)
    return None   

def plot_read_length_pos(
    save_dir:PosixPath,
    load_query:PosixPath,
    load_seq_json:PosixPath,
    smpl:str,
    ref:Ref,
    offset=15,
    plot_region=[]
):
    save_dir = save_dir / 'plot_read_length_pos'
    if not save_dir.exists():
        save_dir.mkdir()
    
    df_query = pd.read_csv(load_query,index_col=0,header=0)
    reporter = df_query.index[0].split('_0')[0]
    dict_query = dict(zip(list(df_query.index),df_query.iloc[:,0].values))
    with gzip.open(load_seq_json,'rt') as f:
        dict_seq = json.load(f)
    
    seq_reporter = ref.seq_reporters[reporter]
    l_reporter = len(seq_reporter)
    
    dict_seq_length = {}
    pos_len_mat = np.zeros((80,l_reporter))
    colors = []
    for name_query, seq_query in dict_query.items():

        l = len(seq_query)
        pos = int(name_query.split('_')[-1])
        if name_query not in dict_seq:
            continue
        length_list = []
        for seq in dict_seq[name_query]:
            # trim GGG and polyA
            if seq.find('GGG') == 0:
                seq = seq[3:]
            polyA_pos = seq.find('AAAAAAA')
            if polyA_pos < 0:
                tmp = find_poly_anywhere(
                        seq=seq,
                        base = 'A',
                        min_len = 7,
                        min_frac = 0.7
                    )
                if tmp is not None:
                    polyA_pos = tmp[0]
                else:
                    polyA_pos = -1
            if polyA_pos > 0:
                seq = seq[:polyA_pos]
                colors.append(['#000000'])
            else:
                print("No polyA found")
                colors.append(['#707070'])
            l_seq = len(seq)
            pos_len_mat[l_seq,pos+offset] += 1
            length_list.append(l_seq)
        dict_seq_length[name_query] = length_list
    
    # save
    with gzip.open(save_dir / f'read_length_pos_{reporter}_{smpl}.json.gz','wt') as f:
        json.dump(dict_seq_length,f)

    idx_nonzero = np.where(np.any(pos_len_mat,axis=1))[0]
    mat = pos_len_mat[idx_nonzero[0]:idx_nonzero[-1]+1,:]
    if len(plot_region) == 2:
        mat = mat[:,plot_region[0]:plot_region[1]]
    else:
        plot_region = [0,l_reporter]

    # plot as heatmap
    fig, ax = plt.subplots(1,1,figsize=(15,3))

    sns.heatmap(
        mat,
        ax=ax,
        cmap='Blues',
        cbar_kws={'label':'Counts'},
        vmin=0,vmax=mat.max()/2
    )
    ax.set_xlabel(f'Position along the reporter (nt), offset={offset} nt')
    ax.set_ylabel('Read length (nt)')
    ax.set_title(smpl, loc='left') 
    ax.set_yticks(np.arange(0,idx_nonzero[-1]-idx_nonzero[0]+1,5))
    ax.set_yticklabels(np.arange(idx_nonzero[0],idx_nonzero[-1]+1,5),rotation=0)
    ax.set_xticks(np.arange(0,mat.shape[1],50))
    ax.set_xticklabels(
        np.arange(0,mat.shape[1],50) + plot_region[0] if len(plot_region) == 2 else np.arange(0,mat.shape[1],50)
        ,rotation=0)

    for region_dict in ref.ref_reporters[reporter]['features']:
        region_name = region_dict['name']
        if 'rich' in region_name:
            color = '#2deade'
        
        if len(plot_region) == 2:
            if region_dict['end'] < plot_region[0]:
                continue
            if region_dict['start'] > plot_region[1]:
                continue
            if region_dict['start'] < plot_region[0]:
                region_dict['start'] = plot_region[0]
            if region_dict['end'] > plot_region[1]:
                region_dict['end'] = plot_region[1]

        ax.axvspan(
            region_dict['start'] - plot_region[0],
            region_dict['end'] - plot_region[0],
            color=region_dict['color'],
            alpha=0.25, lw=0, zorder=10)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(save_dir / f'read_length_pos_{reporter}_{smpl}.pdf')
    plt.close()    

def plot_indv_subsequence(
    save_dir:PosixPath,
    load_query:PosixPath,
    load_seqcnt:PosixPath,
    reporter:str,
    smpl:str,
    ref:Ref=None,
    offset=15,
    ylim=[],
    plot_region=[],
    pause_score_bg = [],
    pause_score_target = [],
    suffix=''
):
    save_dir = save_dir / 'plot_subsequence'
    if not save_dir.exists():
        save_dir.mkdir()
    
    df_query = pd.read_csv(load_query,index_col=0,header=0)
    dict_query = dict(zip(list(df_query.index),df_query.iloc[:,0].values))
    
    df = pd.read_csv(load_seqcnt,header=0,index_col=0)
    df = df[ df.index > 0 ]
    if df.shape[1] > 1:
        df = df.loc[:,[smpl]]
    # df['pos'] = [
    #     ref.seq_reporters[reporter].find(dict_query[f'{reporter}_{pos}']) + offset
    #     for pos in df.index
    # ]
    df['pos'] = np.array(list(df.index)) + offset
    df.reset_index(drop=False,inplace=True)
    df = df.set_index('pos')
    df.drop(['index'],axis=1,inplace=True)
    fname = load_seqcnt.stem.split('.')[0]

    if len(plot_region) == 0:
        if ref is not None:
            plot_region = [0,len(ref.seq_reporters[reporter]) + offset]
        else:
            plot_region = [0,df.index.max()+1]

    avg_counts = df.mean(axis=0)
    df_all = df.copy()
    df = df.loc[plot_region[0]:plot_region[1]]
    
    ### plot
    outfile_name = save_dir / f'indiv_plot{suffix}_{fname}'
    pdf = PdfPages(outfile_name.with_suffix(".pdf")) 

    fig,ax = plt.subplots(1,1,figsize=(4,2))
    # x = df.values / avg_counts.mean()
    # normalize by the first Pro peak
    x_norm = df.loc[pause_score_target[0]-6:pause_score_target[0]-4].values
    if len(x_norm.shape) == 2:
        x_norm = x_norm.mean(axis=1)
    x_norm = np.sum(x_norm)
    if x_norm == 0:
        x_norm = 1
    x = df.values / x_norm
    if len(x.shape) == 2:
        x = x.mean(axis=1)
    ax.plot(
        x,
        color='#808080',marker='.'
    )
    # ax.axhline(1,color='#ff0000',ls='--',lw=0.5)
    seq = ref.seq_reporters[reporter][plot_region[0] : plot_region[1]]
    codon_seq = my._seq2codon(seq)
    if len(ylim) > 0:
        pos = -ylim[0]*0.5
    else:
        pos = np.amax(x) * 0.1
    
    aa_seq = my._seq2aa(seq)
    for j in range(len(aa_seq)):
        ax.text(
            3*j+1,-pos,aa_seq[j],
            ha='center',va='top',fontsize=6,
        )
    if len(ylim) > 0:
        ax.set_ylim(bottom=ylim[0],top=ylim[1])
    else:
        ax.set_ylim(bottom=-pos*3,top=np.amax(x)*1.1)
    
    ax.set_xticks(np.arange(0,plot_region[1]-plot_region[0],20))
    ax.set_xticklabels(
        np.arange(0,plot_region[1]-plot_region[0],20) + plot_region[0],
        rotation=0
    )
    ax.set_xlabel(f'Position along the {reporter} (nt)')
    ax.set_ylabel(f'Read density\n(relative to first Pro)')
    # ax.set_ylabel(f'Read density\n(relative to Pro {pause_score_target[0]-6}-{pause_score_target[0]-3} nt)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(pdf,format="pdf")
    plt.close()

    # plot the full density
    for i in range(2):

        fig,ax = plt.subplots(1,1,figsize=(8,1))
        x = df_all.values / avg_counts.mean()
        if len(x.shape) == 2:
            x = x.mean(axis=1)
        sns.heatmap(
            x.reshape(1,-1),
            ax=ax,
            cmap='Blues',
            cbar=False,
            vmin=0,vmax=10
        )
        # annotate the features
        if i == 0:
            if ref is not None:
                for region_dict in ref.ref_reporters[reporter]['features']:
                    region_name = region_dict['name']
                    if 'rich' in region_name:
                        color = '#2deade'
                    ax.axvspan(
                        region_dict['start'] - df_all.index[0],
                        region_dict['end'] - df_all.index[0],
                        color=region_dict['color'],
                        alpha=0.5, lw=0, zorder=10)
                    
                    if 'start' in region_name.lower():
                        start_pos = region_dict['start']
                    if 'stop' in region_name.lower():
                        stop_pos = region_dict['start']

        ax.set_xlabel(f'Position along the {reporter} (nt)')
        ax.set_xticks(np.arange(0,df_all.index.max()-df_all.index.min(),100))
        ax.set_xticklabels(
            np.arange(0,df_all.index.max()-df_all.index.min(),100),
            rotation=0
        )
        # off y-ticks
        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.tight_layout()
        fig.savefig(pdf,format="pdf")
        plt.close()
    
    pdf.close()

    if len(pause_score_bg) == 2 & len(pause_score_target) == 2:
        # compute pause score using P2A region
        x = df_all.loc[pause_score_target[0]:pause_score_target[1]-1].values.flatten()
        y = df_all.loc[pause_score_bg[0]:pause_score_bg[1]-1].values.flatten()
        frame_x = (np.arange(pause_score_target[0],pause_score_target[1]) - start_pos) % 3
        frame_y = (np.arange(pause_score_bg[0],pause_score_bg[1]) - start_pos) % 3

        ratio = x.sum() / y.sum()
        ratio_frame0 = x[ frame_x==0 ].sum() / y[ frame_y==0 ].sum()

        seq = ref.seq_reporters[ reporter ]
        seq_target = seq[ pause_score_target[0]:pause_score_target[1] ]
        seq_bg = seq[ pause_score_bg[0]:pause_score_bg[1] ]

        with open(save_dir / f'counts{suffix}.txt','w') as f:
            # write sequence + count for target
            f.write(f'# pause score (all frames) x/y = {ratio:.3f}\n')
            f.write(f'# pause score (frame 0) x/y = {ratio_frame0:.3f}\n')
            f.write(f'# x ({pause_score_target[0]}:{pause_score_target[1]} nt)\n')
            f.write(' sequence\tcount\tframe\n')
            for i in range(len(x)):
                f.write(f'{seq_target[i]}\t{int(x[i])}\t{int(frame_x[i])}\n')
            f.write(f'# y ({pause_score_bg[0]}:{pause_score_bg[1]} nt)\n')
            f.write(' sequence\tcount\tframe\n')
            for i in range(len(y)):
                f.write(f'{seq_bg[i]}\t{int(y[i])}\t{int(frame_y[i])}\n')
            
            






def pause_score(
    save_dir:PosixPath,
    dict_load_seqcnt:dict,
    dict_columns:dict,
    dict_reporters:dict,
    dict_pos_target:dict,
    dict_pos_bg:dict,
    ref:Ref=None,
    offset=15,
    width=3,height=3,
    suffix=''
):
    save_dir = save_dir / 'plot_subsequence'
    if not save_dir.exists():
        save_dir.mkdir()
    
    seq_P2A = 'gctactaacttcagcctgctgaagcaggctggagacgtggaggagaaccctggacct'.upper()
    
    pause_scores = {}
    for smpl, load_seqcnt in dict_load_seqcnt.items():
        df = pd.read_csv(load_seqcnt,header=0,index_col=0)
        df = df.loc[:,[dict_columns[smpl]]]
        df['pos'] = np.array(list(df.index)) + offset
        df.reset_index(drop=False,inplace=True)
        df = df.set_index('pos')
        df.drop(['index'],axis=1,inplace=True)
        x = df.loc[dict_pos_target[smpl][0]:dict_pos_target[smpl][1]].sum().iloc[0]
        y = df.loc[dict_pos_bg[smpl][0]:dict_pos_bg[smpl][1]].sum().iloc[0]

        seq = ref.seq_reporters[ dict_reporters[smpl] ]
        seq_target = seq[ dict_pos_target[smpl][0] : dict_pos_target[smpl][1] ]
        seq_bg = seq[ dict_pos_bg[smpl][0] : dict_pos_bg[smpl][1] ]

        assert seq_bg.upper() == seq_P2A

        pause_scores[smpl] = x / y
    
    # bar plot
    outfile_name = save_dir / f'pause_scores{suffix}'
    pdf = PdfPages(outfile_name.with_suffix(".pdf")) 

    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.bar(
        pause_scores.keys(),
        pause_scores.values(),
        color='#353535'
    )
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(pause_scores.keys(),rotation=45,ha='right')
    ax.set_ylabel('Pause score')
    ax.set_ylim(0,0.7)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(pdf,format="pdf")
    plt.close()

    pdf.close()

        
