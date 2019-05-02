'''
from Singhal et al 2015

https://github.com/singhal/postdoc/blob/master/simulations_hotspot_power.py
'''

import numpy as np
import re
import sys
from itertools import zip_longest as izip
import subprocess
import random
import os
import argparse
from tqdm import tqdm
from datetime import datetime

def args():
    parser = argparse.ArgumentParser(description = 'run macs and process output',
                                     usage = 'python3.5 simulations_hotspot_power.py [options]')

    parser.add_argument('-o', '--outdir', required = True,
                        type = str, help = 'Dir to write output to')

    args = parser.parse_args()

    return args.outdir


def simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, diffs, hotspot_lengths, ix):
    
    print('writing hotspot file')
    hotspot_file = '%shotspots/hotspot_rho%s_%s.txt' % (out_dir, rho, ix)
    hotspot_f = open(hotspot_file, 'w')
    num_hotspots = len(diffs) * len(hotspot_lengths)
    starts = np.linspace(0+50e3,seq_size-50e3,num_hotspots)
    for i, diff in enumerate(diffs):
        for j, length in enumerate(hotspot_lengths):
            spot_start = starts[j + i * len(hotspot_lengths)]
            spot_end = spot_start + length
            hotspot_f.write('%.4f\t%.4f\t%.5f\n' % (spot_start / float(seq_size), spot_end / float(seq_size), diff))
    hotspot_f.close()    

    haplo_file = '%shaplotypes/haplo_rho%s_%s.fa' % (out_dir, rho, ix)
    anc_file = '%sancallele/ancallele_rho%s_%s.txt' % (out_dir, rho, ix)
    haplo_f = open(haplo_file, 'w')
    anc_f = open(anc_file, 'w')

    with open('macs.err', 'w') as f_err:
        macs = subprocess.Popen('./bin/macs %s %s -t %s -r %s -R %s | bin/msformatter' % (nsam, seq_size, theta, rho, hotspot_file), shell=True, stdout=subprocess.PIPE, stderr=f_err)

    positions = []
    haplo = []
    for l in macs.stdout:
        l = l.decode('utf-8')
        if re.match(r'^positions:', str(l)):
            positions = [float(match) for match in re.findall('([\d\.e-]+)', l)]
        haplo.append(l.rstrip())            

    haplo = haplo[(len(haplo) - nsam):len(haplo)]
    for ix, hap in enumerate(haplo):
        haplo[ix] = list(hap)
    # these don't have singletons
    singletons = []
    positions_pruned = []
    haplo_pruned = []
    print('finding singletons')
    for ix, pos in tqdm(enumerate(positions)):
        allele_count = 0
        for hap in haplo:
            if hap[ix] == '1':
                allele_count += 1
        if allele_count <= 1:
            singletons.append(ix)

    for ix, pos in tqdm(enumerate(positions)):
        if ix not in singletons:
            pos = int(round(pos * seq_size)) - 1 
            if pos < 0:
                pos = 0
            if pos > (seq_size -1):
                pos = seq_size - 1
            positions_pruned.append(pos)
    
    print('removing singletons')
    for ix1, hap in tqdm(enumerate(haplo)):
        haplo_pruned.append([])
        for ix2, pos in enumerate(hap):
            if ix2 not in singletons:
                haplo_pruned[ix1].append(pos)

    del positions
    del haplo
    del singletons
    
    print('assigning bases')
    bases = []
    for base, freq in eq_freq.items():
        bases += [base] * int(freq * 1000)
    seq = []
    for x in range(int(seq_size)):
        seq.append(random.choice(bases))

    mutations = {}
    
    print('assigning muts')
    for pos in positions_pruned:
        anc = seq[pos]
        ran_num = random.random()
        for base, prob in mut_rates[anc].items():
            ran_num = ran_num - prob
            if ran_num < 0:
                mut = base
                nuc = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
                prob = [0.02, 0.02, 0.02, 0.02]
                prob[nuc[mut]] = 0.94
                anc_f.write('%s %s\n' % (pos, ' '.join('%.2f' % i for i in prob)))
                mutations[pos] = [anc, mut]
                break
    anc_f.close()

    print('writing haplotypes...')
    for ind, ind_hap in tqdm(enumerate(haplo_pruned)):
        tmp_seq = seq[:]
        for hap_ix, bp in enumerate(ind_hap):
            if bp == '1':
                tmp_seq[positions_pruned[hap_ix]] = mutations[positions_pruned[hap_ix]][1]
        haplo_f.write('>haplo%s\n' % ind)
        for seq_i in range(0, len(tmp_seq), 60):
            haplo_f.write('%s\n' % ''.join(tmp_seq[seq_i:seq_i+60]))
    haplo_f.close()


def main():
    out_dir = args()
    seq_size = 1000000
    num_sim = 10
    # hotspot / coldspot difference, > 1
    diffs = [10, 10, 20, 20, 40, 40, 60, 60]
    # mean rho values
    rhos = [0.0001, 0.001, 0.01, 0.1, 1.0, 2.5]
    # rhos = [0.0001]
    # hotspot length
    hotspot_lengths = [2000]
    # theta (per bp!)
    theta = 0.03
    # number of haplotypes to sample
    nsam = 19
    # A, C, T, G
    eq_freq = {'A': 0.179, 'C': 0.321, 'G': 0.320, 'T': 0.180}
    # mutation rates, modified from mutation matrix
    mut_rates = {'A': {'C': 0.282, 'G': 0.487, 'T': 0.231},
                 'C': {'A': 0.260, 'G': 0.220, 'T': 0.520},
                 'G': {'A': 0.515, 'C': 0.227, 'T': 0.258},
                 'T': {'A': 0.262, 'C': 0.476, 'G': 0.262}}
    
    for rho in rhos:
        for ix in range(100, 100+num_sim):
            print('currently on {rho}, sim {ix}'.format(rho=rho, ix=ix))
            haplo_file = '%shaplotypes/haplo_rho%s_%s.fa' % (out_dir, rho, ix)
            if not os.path.isfile(haplo_file):
                simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, diffs, hotspot_lengths, ix)
                print('completed {rho}, sim {ix}'.format(rho=rho, ix=ix))
                print('currently', datetime.now())
    

if __name__ == "__main__":
    main()
