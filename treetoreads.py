#!/usr/bin/env python
"""Tree to Reads - A python script to to read a tree,
    resolve polytomies, generate mutations and simulate reads."""
import random
import os
import sys
import argparse
from subprocess import call
import dendropy
import itertools
from operator import itemgetter

VERSION = "0.0.5"

tmpout = open("tmpdels.txt","w")

class TreeToReads(object):
    """A tree to reads object that holds the input tree and base genome,
    and has methods to create different outputs"""
    _argread = 0
    _treeread = 0
    _simran = 0
    _madeout = 0
    _siteread = 0
    _mutlocs = 0
    _genmut = 0
    _genread = 0
    _vargen = 0
    def __init__(self, configfi, run=1, main=None):
        """initialized object, most attributes generated through self._check_args using config file."""
        #self.seed = (5)
        self.seed = random.randint(0, sys.maxint)
       # self.seed = (4871050314122165556)
        sys.stdout.write("Random seed is {}\n".format(self.seed))
        random.seed(self.seed)
        self.configfi = configfi
        self.run = run
        self.main = main
        if os.path.isfile(self.configfi):
            sys.stdout.write("Running TreetoReads using configuration file {}\n".format(self.configfi))
        else:
            sys.stderr.write("Config file '{}' not found. Exiting.\n".format(self.configfi))
            self._exit_handler()
        self._check_args()
        self.test_deps()
        if 'prefix' in self.config:
            self.prefix = self.config['prefix']
        else:
            self.prefix = ''
        if self.run:
            if (self.no_art == 0) and (self.config.get('coverage')):
                self.run_art()
            else:
                sys.stdout.write("simulating genomes but not reads\n")
                if self.get_arg("indel_model"):
                    if call(['which', 'indelible'], stdout=open('/dev/null', 'w')) == 1:
                        sys.stderr.write('''ERROR: indelible not found. Needs to be installed to simulate insertions
                                         and deletions. IGNORING indel parameters and continuing simulation\n''')
                        self.mut_genomes_no_indels()
                    else:
                        self.mut_genomes_indels()
                else:
                    self.mut_genomes_no_indels()
    def _exit_handler(self):
        '''makes debugging interactively easier, by not exiting on errors'''
        if self.main:
            sys.exit()
        else:
            pass
    def test_deps(self):
        "Check that seq-gen and ART are installed,"
        if call(['which', 'seq-gen'], stdout=open('/dev/null', 'w')) == 1:
            sys.stderr.write('''seq-gen needs to be installed
                              and in your path for TreeToReads to run.
                              It was not found.  Exiting\n''')
            self._exit_handler()
        if call(['which', 'art_illumina'], stdout=open('/dev/null', 'w')) == 1:
            sys.stderr.write('''ERROR: art_illumina needs to be installed
                             and in your path for TreeToReads to
                             generate reads. Art not found.
                             TTR will only generate mutated genomes. \n''')
            self.no_art = 1
        else:
            self.no_art = 0
        if not dendropy.__version__.startswith('4'):
            sys.stderr.write('''ERROR: Please upgrade the python package dendropy to version 4,
                                using 'pip install dendropy --upgrade'.
                                Exiting\n''')
            self._exit_handler()

    def read_args(self):
        """reads arguments from config file"""
        try:
            config = open(self.configfi)
        except:
            sys.stderr.write('''Config file '{}' not found.
                              Exiting.'''.format(self.configfi))
            self._exit_handler()
        self.config = {}
        poss_args = ['treefile_path',
                     'number_of_variable_sites',
                     'base_genome_name',
                     'base_genome_path',
                     'rate_matrix',
                     'freq_matrix',
                     'coverage',
                     'prefix',
                     'output_dir',
                     'mutation_clustering',
                     'percent_clustered',
                     'exponential_mean',
                     'gamma_shape',
                     'read_length',
                     'fragment_size',
                     'stdev_frag_size',
                     'error_model1',
                     'error_model2',
                     'indel_model',
                     'indel_rate']
        for lin in config:
            lii = lin.split('=')
            self.config[lii[0].strip()] = lii[-1].split('#')[0].strip()
            if (lii[0].strip() != '') and (not lii[0].strip().startswith("#")) and (lii[0].strip() not in poss_args):
                sys.stderr.write('''Config paramater '{}' not in possible parameter list. Acceptable params are {}
                                 Exiting\n.'''.format(lii[0].strip(), "\n".join(poss_args)))
                self._exit_handler()
        sys.stdout.write("Arguments read\n")
        self._argread = 1

    def make_output(self):
        """Creates output directory"""
        if not self._argread:
            self.read_args()
        self._madeout = 1
        sys.stdout.write('output directory is {}\n'.format(self.outd))
        if not os.path.isdir(self.outd):
            os.mkdir(self.outd)
        configout = open('{}/analysis_configuration.cfg'.format(self.outd), 'w')
        for lin in open(self.configfi).readlines():
            configout.write(lin)
        configout.close()

    def get_arg(self, nam):
        """Returns arugments from the argument dictionary"""# a bit more convoluted than necessary...
        if not self._argread:
            self.read_args()
        if nam in self.argdict:
            return self.config[self.argdict[nam]]
        else:
            return self.config.get(nam)

    def _check_args(self):
        """Checks that arguments are of the appropriate types,
        and all required args are present."""
        if self._argread != 1:
            self.read_args()
        self.argdict = {'treepath':'treefile_path',
                        'nsnp':'number_of_variable_sites',
                        'base_name':'base_genome_name',
                        'genome':'base_genome_path',
                        'ratemat':'rate_matrix',
                        'cov':'coverage',
                        'outd':'output_dir'
                       }
        for arg in self.argdict:
            if self.argdict[arg] not in self.config:
                sys.stderr.write("{} is missing from the config file".format(self.argdict[arg]))
                self._exit_handler()
        try:
            self.nsnp = int(self.get_arg('nsnp'))
            sys.stdout.write('Number of variable sites is {}\n'.format(self.nsnp))
        except:
            sys.stderr.write('''number of variable {} could not be coerced to an integer.
                              Exiting.\n'''.format(self.get_arg('nsnp')))
            self._exit_handler()
        if len(self.get_arg('ratemat').split(',')) == 6:
            self.ratemat = dict(zip(['ac', 'ag', 'at', 'cg', 'ct', 'gt'],self.get_arg('ratemat').split(',')))
        else:
            sys.stderr.write('''{} values in rate matrix, there should be 6.
                              Exiting.\n'''.format(len(self.get_arg('ratemat').split(','))))
            self._exit_handler()
        try:
            tf = open(self.get_arg('treepath'))
            if tf.readline().startswith('#NEXUS'):
                self.treetype = "nexus"
            else:
                self.treetype = "newick"
            tf.close()
        except:
            sys.stderr.write('''Could not open treefile {}.
                              Exiting.\n'''.format(self.get_arg('treepath')))
            self._exit_handler()
        try:
            open(self.get_arg('genome'))
        except:
            sys.stderr.write('''Could not open base genome {}.
                              Exiting.\n'''.format(self.get_arg('genome')))
            self._exit_handler()
        try:
            self.outd = self.get_arg('outd')
        except:
            self.outd = ('ttr_out')
            sys.stderr.write("Setting output directory to ttr_out\n")
        #optional parameters
        if self.get_arg('gamma_shape') is not None:
            try:
                self.shape = float(self.get_arg('gamma_shape'))
            except:
                sys.stderr.write('''shape parameter {}
                                 could not be coerced to a float.
                                 Exiting.\n'''.format(self.get_arg('gamma_shape')))
                self._exit_handler()
        else:
            self.shape = None
        if self.get_arg('mutation_clustering') is not None:
            if self.get_arg('mutation_clustering') == 'ON':
                self.clustering = 1
                try:
                    self.clust_perc = float(self.get_arg('percent_clustered'))
                    self.lambd = 1.0/(int(self.get_arg('exponential_mean'))-1)
                    sys.stdout.write("clustering proportion is {}\n".format(self.clust_perc))
                    sys.stdout.write("exponential_mean is {}\n".format(self.get_arg('exponential_mean')))
                except:
                    sys.stderr.write('''Problem reading clustering parameters,
                                     requires float for 'percent_clustered'
                                     and 'exponential_lambda'. Exiting.\n''')
                    self._exit_handler()
            else:
                sys.stdout.write('''Mutation clustering is OFF,
                                  to use set mutation_clustering = ON 
                                  and values for "percent_clustered" 
                                  and an integer for "exponential_mean"\n''')
                self.clustering = 0
        else:
            sys.stdout.write('Mutation clustering is OFF\n')
            self.clustering = 0

    def read_tree(self):
        """Reads in a tree from a file, arbitrarily resolves poltomies if present,
        strips leading [&U] and writes out to outputdir/simtree.tre"""
        self._treeread = 1
        if not self._madeout:
            self.make_output()
        if not self._genread:
            self.read_genome
        #import tree from path
        if dendropy.__version__.startswith('4'):
            taxa = dendropy.TaxonNamespace()
            try:
                tree = dendropy.Tree.get_from_path(self.get_arg('treepath'), self.treetype, taxon_namespace=taxa, preserve_underscores=True)
            except:
                sys.stderr.write("Problems reading the tree - is it in proper newick or nexus format?\n")
                self._exit_handler()
        else:
            taxa = dendropy.TaxonSet()
            try:
                tree = dendropy.Tree.get_from_path(self.get_arg('treepath'), self.treetype, taxon_set=taxa, preserve_underscores=True)
            except:
                sys.stderr.write("Problems reading the tree - is it in proper newick or nexus format?\n")
                self._exit_handler()
        if tree.length() == 0:
            sys.stderr.write("TTR requires branch lengths. Branch lengths appear to be missing (treelength = 0). Exiting.\n")
            self._exit_handler()
        self.seqnames = taxa.labels()
        self.base_name = self.get_arg('base_name')
        if self.base_name not in self.seqnames:
            sys.stderr.write("base genome name {} is not in tree. Exiting.\n".format(self.base_name))
            self._exit_handler()
        tree.resolve_polytomies()
        tree_len = tree.length()
        expected_tree_len = float(self.nsnp)/self.genlen
        for edge in tree.postorder_edge_iter():
            if edge.length is None:
                edge.length = 0
            else:
                edge.length = (float(edge.length)/tree_len) * expected_tree_len
        assert -0.001 < expected_tree_len - tree.length() < 0.001
        self.scaledouttree = "{}/scaledtree.tre".format(self.outd)
        tree.write_to_path(self.scaledouttree,
                           schema='newick',
                           suppress_internal_node_labels=True,
                           suppress_rooting=True)
        self.scaled_tree_newick = tree.as_string(schema='newick', real_value_format_specifier='.15f')
        if expected_tree_len < 0.01:  #scale up tree length so generate mutations in seqgen without a million invariant sites.
            stretch = 0.01/expected_tree_len
            for edge in tree.postorder_edge_iter():
                if edge.length is None:
                    edge.length = 0
                else:
                    edge.length = edge.length * stretch
        self.outtree = "{}/simtree.tre".format(self.outd)
        tree.write_to_path(self.outtree, schema='newick', suppress_internal_node_labels=True, suppress_rooting=True)
        sys.stdout.write("Tree read\n")

    def read_genome(self):
        """Reads in base genome information to use for simulations from file"""
        self._genread = 1
        if not self._argread:
            self.read_args()
        self.genlen = 0
        contigs = 0
        self.contig_breaks = []
        self.contig_names = []
        base_counts = {'A':0, 'C':0, 'G':0, 'T':0}
        with open(self.get_arg('genome'), 'r') as in_file:
            for line in in_file:
                line = line.strip()
                if line.startswith('>'):
                    contigs += 1
                    self.contig_breaks.append(self.genlen)
                    self.contig_names.append(line.strip('>'))
                else:
                    self.genlen += len(line)
                    if not set(line.upper()).issubset(set(['A', 'T', 'G', 'C', 'N'])):
                        sys.stderr.write('''Your genome appears to have characters other than ATGC,
                                         such as: {} Please check your input genome.\n'''.format(set(line)))
                    for base in base_counts:
                        base_counts[base] += line.count(base)
        self.contig_breaks.append(self.genlen)
        self.freqmat = {base:str(base_counts[base]/float(self.genlen)) for base in ['A','C','G','T']}
        sys.stdout.write("Base frequencies detected from anchor genome A:{} C:{} G:{} T:{}\n".format(self.freqmat['A'],
                                                                                             self.freqmat['C'],
                                                                                             self.freqmat['G'],
                                                                                             self.freqmat['T']))
        if self.nsnp > self.genlen:
            sys.stderr.write('''number of variables sites {}
                             is higher than the length
                             of the contig or geonme {}.
                             Exiting\n'''.format(self.nsnp, self.genlen))
            self._exit_handler()
        sys.stdout.write("genome has {} contigs\n".format(len(self.contig_breaks)))
        sys.stdout.write("Genome has {} bases\n".format(self.genlen))

    def generate_varsites(self):
        """Runs seqgen to generate variable sites on tree
        seqgen GTR specified as : A to C, A to G, A to T, C to G, C to T and G to T
        and freqs, seqgen uses a,c,g,t UNLIKE INdelible
        """
        self._vargen = 1
        if not self._treeread:
            self.read_tree()
        self.simloc = "{}/seqs_sim.txt".format(self.outd)
        lenseqgen = 120*self.nsnp
        freqs = ','.join([self.freqmat[base] for base in ['A','C','G','T']])
        gtr = ','.join([self.ratemat[rate] for rate in ['ac','ag','at','cg','ct','gt']])
        if self.shape:
            seqgenpar = ['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR',
                         '-a{}'.format(self.shape), '-r{}'.format(gtr),
                         '-f{}'.format(freqs), '-or']
        else:
            seqgenpar = ['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR',
                         '-r{}'.format(gtr),
                         '-f{}'.format(freqs), '-or']
        call(seqgenpar,
             stdout=open('{}'.format(self.simloc), 'w'),
             stderr=open('{}/seqgen_log'.format(self.outd), 'w'),
             stdin=open('{}'.format(self.outtree)))
        assert  open('{}/seqgen_log'.format(self.outd)).readlines()[-1].startswith("Time taken")
        sys.stdout.write("Variable sites generated using seq-gen\n")

    def read_varsites(self, add=False):
        """Reads in only the variable sites from the seqgen output file
        Stores as ditionary."""
        self._siteread = 1
        if not self._vargen:
            self.generate_varsites()
        if add is False:
            self.sitepatts = {}
            for nuc in ['A', 'G', 'T', 'C']:
                self.sitepatts[nuc] = []
            self.trip_hit = 0
            self.var_site = 0
        nucsets = {}
        varsites = []
        with open(self.simloc) as f:
            next(f)
            for lin in f:
                bases = lin.split()[1].strip()
                for i, nuc in enumerate(bases):
                    if i not in nucsets:
                        nucsets[i] = set()
                    nucsets[i].add(nuc)
        for i in nucsets:
            if len(nucsets[i]) >= 2:
                varsites.append(i)
                self.var_site += 1
            if len(nucsets[i]) > 2:
                self.trip_hit += 1
        simseqs = {}
        with open(self.simloc) as f:
            next(f)
            for lin in f:
                seq = lin.split()[0]
                if seq.startswith('"') and seq.endswith('"'):
                    seq = seq[1:-1]
                if seq.startswith("'") and seq.endswith("'"):
                    seq = seq[1:-1]
                simseqs[seq] = []
                bases = lin.split()[1].strip()
                for i in varsites:
                    simseqs[seq].append(bases[i])
        try:
            assert set(self.seqnames) == set(simseqs.keys())
        except:
            sys.stderr.write('''Seqnames don't match simulated keys,
                             seqnames :{} simulated names: {}
                             '''.format('\n'.join(self.seqnames), '\n'.join(simseqs.keys())))
            self._exit_handler()
        ref = simseqs[self.base_name]
        for i, nuc in enumerate(ref):
            site = {}
            nucs = set()
            for srr in simseqs:
                site[srr] = simseqs[srr][i]
                nucs.add(simseqs[srr][i])
            assert len(nucs) > 1
            self.sitepatts[nuc].append(site)

    def add_varsites(self):
        """Simulates more varible sites if there are not enough"""
        sys.stderr.write('simulating additional variable sites\n')
        if self._siteread != 1:
            self.read_varsites()
        self.generate_varsites()
        self.read_varsites(add=True)


    def select_mutsites(self):
        """Selects which positions in the base genome will be variable sites."""
        if not self._madeout:
            self.make_output()
        if not self._genread:
            self.read_genome()
        self.mutsite = "{}/mutsites.txt".format(self.outd)
        fi = open(self.mutsite, "w")
        rands = set()
        if self.clustering:
            nclust = int((self.nsnp*self.clust_perc)/2)
            ranpairA = random.sample(range(self.genlen), nclust)
            for site in ranpairA:
                rands.add(site)
                diff = int(random.expovariate(self.lambd)) + 1
                if (random.choice([0, 1]) or (site-diff < 0)) and (site+diff < self.genlen):
                    ranpairB = site+diff
                else:
                    ranpairB = site-diff
                rands.add(ranpairB)
            ransingle = random.sample(range(self.genlen), int(self.nsnp*(1-self.clust_perc)))
            rands.update(ransingle)
        else:
            ran = random.sample(range(self.genlen), self.nsnp)
            rands = set(ran)
        while len(rands) < self.nsnp: #deals inelegantly with multiple hits, to make sure there are nsnp-len individual sites
            ran = random.sample(range(self.genlen), (self.nsnp-len(rands)))
            rands = rands | set(ran)
        self.mutlocs = list(rands)
        self.mutlocs.sort()
        

    def assign_sites(self):
        """Pulls columns for each mutation"""
        self.mutations = {}
        if not self._mutlocs:
            self.select_mutsites()
        if not self._siteread:
            self.read_varsites()
        self.mut_genos = {}
        patnuc = {}
        ii = 0
        patnuc['A'] = 0
        patnuc['G'] = 0
        patnuc['T'] = 0
        patnuc['C'] = 0
        self.snpdic = {}
        #This next section is to index mutations, and re-simulate if there are not enough of a base
        with open(self.get_arg('genome'), 'r') as in_file:
            for line in in_file:
                if line.startswith('>'):
                    pass
                else:
                    line = line.strip()
                    for nuc in line:
                        if ii in self.mutlocs:
                            patnuc[nuc] += 1
                            self.snpdic[ii] = patnuc[nuc]
                            if len(self.sitepatts[nuc]) <= self.snpdic[ii]:
                                self.add_varsites()
                        ii += 1
        if self.trip_hit:
            sys.stderr.write('''{} variable sites with more than 2 bases out of {} sites.
            Scale down tree length if this is too high.\n'''.format(self.trip_hit, self.var_site))


    def mut_genomes_no_indels(self):
        """Writes out the simulated genomes with mutations"""
        self.assign_sites()
        matout = open("{}/var_site_matrix".format(self.outd), 'w')
        self.vcf_dict = {}
        for loc in self.mutlocs:
            self.vcf_dict[loc] = {}
        self.mut_trans = {}
        contig = 0
        for loc in self.mutlocs:
            adjusted_loc = loc - self.contig_breaks[contig]
            if self.contig_breaks[contig+1] < loc:
                contig += 1
            adjusted_loc = loc - self.contig_breaks[contig]
            contig_name = self.contig_names[contig]
            self.mut_trans[loc] = (contig_name, str(adjusted_loc))
        for seq in self.seqnames:
            self.mut_genos[seq] = []
            sys.stdout.write("writing genome for {}\n".format(seq))
            if not os.path.isdir("{}/fasta_files".format(self.outd)):
                os.mkdir("{}/fasta_files".format(self.outd))
            genout = open("{}/fasta_files/{}{}.fasta".format(self.outd, self.prefix, seq), 'w')
            ii = 0
            lw = 0 #tracking line wrapping for fasta
            with open(self.get_arg('genome'), 'r') as in_file:
                for line in in_file:
                    if line.startswith('>'):
                        if ii > 0:
                            genout.write('\n')
                        #    genout.write(line.strip()+"_"+self.prefix+seq+"\n")
                            lw = 0
                        #else:
                        genout.write(line.strip()+"_"+self.prefix+seq)
                    else:
                        line = line.strip()
                        for nuc in line:
                            if lw%70 == 0:
                                genout.write('\n')
                            if ii in self.mutlocs:
                                contig_name, adjusted_loc = self.mut_trans[ii]
                                if nuc == 'N':
                                    genout.write('N')
                                    self.vcf_dict[ii][seq] = 'N'
                                else:
                                    patt = self.sitepatts[nuc][self.snpdic[ii]]
                                    genout.write(patt[seq])
                                    self.mut_genos[seq].append(patt[seq])
                                    matout.write("{} {} {} {} {}\n".format(seq, patt[seq], ii, contig_name, adjusted_loc))
                                    self.vcf_dict[ii][seq] = patt[seq]
                            else:
                                genout.write(nuc)
                            ii += 1
                            lw += 1
            genout.write('\n')
            genout.write('\n')
            genout.close()
        matout.close()
        self._genmut = 1
        write_vcf(self)
        sys.stdout.write("Mutated genomes\n")

    def mut_genomes_indels(self):#TODO does not account for SNPs in insertions
        """Writes out the simulated genomes with mutations and indels
         indelible GTR is specified as TC TA TG CT CA CG, CG = 1
         and indelible base frequencies are specified as 
         pi_T, pi_C, pi_A, pi_G"""
        self.assign_sites()
        write_indelible_controlfile(self.outd,
                                    self.ratemat,
                                    self.freqmat,
                                    self.get_arg('indel_model'),
                                    self.get_arg('indel_rate'),
                                    self.scaled_tree_newick[:-20],
                                    self.genlen,
                                    self.seed)
        run_indelible(self.outd)
        self.insertions, self.deletions, self.insertionlocs, self.deletionlist = read_indelible_aln(self)
        matout = open("{}/var_site_matrix".format(self.outd), 'w')
        self.vcf_dict = {}
        for loc in self.mutlocs:
            self.vcf_dict[loc] = {}
        for loc in self.insertionlocs:
            self.vcf_dict[loc] = {}
         #   print loc
       # for seq in self.insertions:
       #     print seq
       #     print self.insertions[seq].keys()
        deletion_ids={}
        self.del_location_dict = {}
        self.vcf_dict_deletions = {}
        for seqname in self.deletions:
            self.del_location_dict[seqname] = {}
            for deletion in self.deletions[seqname]:
                for i in deletion:
                    self.del_location_dict[seqname][i] = deletion
               # if not self.vcf_dict_deletions.get(deletion, None):
               #     self.vcf_dict_deletions[deletion] = {}
               #     self.vcf_dict_deletions[deletion][self.base_name] = ""
               # self.vcf_dict_deletions[deletion][seqname] = ""
     #   for i, loc_tuple in enumerate(self.deletionlist):
     #       self.vcf_dict[loc_tuple[0]] = {}
     #       deletion_ids[loc_tuple] = i
#            del_starts.add(loc_tuple[0])
      #  translate_deletions = {}
        crossmap_dict = {}
        for seq in self.seqnames:
            tmpout.write(seq+"\n")
            tmpout.write(", ".join([str(i) for i in self.deletions[seq]]))
            self.mut_genos[seq] = []
            sys.stdout.write("writing genome for {}\n".format(seq))
            if not os.path.isdir("{}/fasta_files".format(self.outd)):
                os.mkdir("{}/fasta_files".format(self.outd))
            genout = open("{}/fasta_files/{}{}_indel.fasta".format(self.outd, self.prefix, seq), 'w')
            ii = 0 #indexing along reference genome
            ali = 0 #indexing along alignement
            lw = 0 #counting for fasta line wrapping
            with open(self.get_arg('genome'), 'r') as in_file:
                prevnuc = '-'
                for line in in_file:
                    if line.startswith('>'):
                        if ii > 0:
                            genout.write('\n')
                            lw = 0
                        genout.write(line.strip()+"_"+self.prefix+seq)
                    else:
                        line = line.strip()
                        for nuc in line:
                            if lw%70 == 0:
                                genout.write('\n')
                            counted = 0
                            if ii in self.insertionlocs:
                               # #print "insertion happed at {}".format(ii)
                                for pos in self.insertionlocs[ii]:
                                    if seq == self.base_name:
                                        genout.write('-')
                                        self.vcf_dict[ii][seq] = nuc
                                    else:
                                        genout.write(self.insertions[seq][pos])
                                        if self.vcf_dict[ii].get(seq):
                                            self.vcf_dict[ii][seq] += self.insertions[seq][pos]
                                        else:
                                            self.vcf_dict[ii][seq] = nuc + self.insertions[seq][pos]
                                    ali += 1
                                    lw += 1
                                    if lw%70 == 0:
                                        genout.write('\n')
                            if ali in self.del_location_dict[seq]: #This should be exclusive of the columns considered "insertions".
                                genout.write('-')
                                current_del_tuple =self.del_location_dict[seq][ali]
                                if not self.vcf_dict_deletions.get(current_del_tuple, None):
                                    self.vcf_dict_deletions[current_del_tuple] = {} #then create a vfc dict entry for it
                                    for subseq in self.seqnames:
                                        self.vcf_dict_deletions[current_del_tuple][subseq] = prevnuc # and start that entry with the anchor base for each sequence
                                if current_del_tuple not in crossmap_dict: #This should just be the startsite
                                    crossmap_dict[current_del_tuple] = {"alignemntstart":ali, "refloc":ii}
                                    #map the location of the dleetion in the alignemnt to the base in the anchor genome.
                                self.vcf_dict_deletions[current_del_tuple][self.base_name] += nuc
                                self.vcf_dict_deletions[current_del_tuple][seq] += 'x'
                                if ii in self.mutlocs:
                                    self.vcf_dict[ii][seq] = '*' #Deletion at a SNP site
                                ali += 1
                                lw += 1
                                ii += 1
                            elif ii in self.mutlocs:
                                if nuc == 'N':
                                    genout.write('N')
                                    self.vcf_dict[ii][seq] = 'N'
                                else:
                                    patt = self.sitepatts[nuc][self.snpdic[ii]]
                                    genout.write(patt[seq])
                                    self.mut_genos[seq].append(patt[seq])
                                    matout.write("{} {} {}\n".format(seq, patt[seq], ii))
                                    self.vcf_dict[ii][seq] = patt[seq]
                                if not counted:
                                    ali += 1
                                    lw += 1
                                    ii += 1
                            else:
                                if not counted:
                                    genout.write(nuc)
                                    ali += 1
                                    lw += 1
                                    ii += 1
                            prevnuc = nuc
            genout.write('\n')
            genout.write('\n')
            genout.close()
        matout.close()
        #print self.vcf_dict_deletions
        #print crossmap_dict
        overlapping_deletions = {}
        last_dele = (0,)
        for dele in self.deletionlist:
            if (set(dele).intersection(last_dele)):
                if len(set(dele).union(last_dele)) == len(last_dele): #second del is subset
                    seq = self.vcf_dict_deletions[last_dele][self.base_name]
                    dele_pos =  [i for i,x in enumerate(last_dele) if x in set(dele).intersection(last_dele)]
                    dele_seq = ''.join([seq[i] for i in xrange(len(seq)) if i not in dele_pos])
                    assert len(dele_seq) < len(seq)
                    for seqname in self.seqnames:
                        if dele in self.deletions[seqname]:
                            self.vcf_dict_deletions[last_dele][seqname] = dele_seq
                    del self.vcf_dict_deletions[dele] # This means cannot currently cope with 3 overlapping deletions
                else: #second del is subset
                    sys.stderr.write('''Too many overlapping deletions for sensical VCF\n''')
                    exit()
            last_dele = dele
        for dele in self.vcf_dict_deletions:
            refseq = self.vcf_dict_deletions[dele][self.base_name]
            loc = crossmap_dict[dele]['refloc']
            self.vcf_dict[loc] = {}
            for seqn in self.seqnames:
                if seqn == self.base_name:
                    self.vcf_dict[loc][seqn] = self.vcf_dict_deletions[dele][seqn]
                else:
                    if len(self.vcf_dict_deletions[dele][seqn]) == 1:
                        self.vcf_dict[loc][seqn] = refseq
                    else:
                        self.vcf_dict[loc][seqn] = self.vcf_dict_deletions[dele][seqn].rstrip('x')                          
        for seq in self.seqnames: #remove the gaps for later sim.
            out_file = open("{}/fasta_files/{}{}.fasta".format(self.outd, self.prefix, seq), 'w')
            inp = "{}/fasta_files/{}{}_indel.fasta".format(self.outd, self.prefix, seq)
            call(['sed', 's/-//g', inp], stdout=out_file)
        self._genmut = 1
        self.mut_trans = {}
        contig = 0
        for loc in self.vcf_dict.keys():
            adjusted_loc = loc - self.contig_breaks[contig]
            if self.contig_breaks[contig+1] < loc:
                contig += 1
            adjusted_loc = loc - self.contig_breaks[contig]
            contig_name = self.contig_names[contig]
            self.mut_trans[loc] = (contig_name, str(adjusted_loc))
        write_vcf(self)
        sys.stdout.write("Mutated genomes\n")


    def run_art(self, coverage=None):
        """Runs ART to simulate reads from the simulated genomes"""
        if not self._genmut:
            if self.get_arg("indel_model"):
                if call(['which', 'indelible'], stdout=open('/dev/null', 'w')) == 1:
                    sys.stderr.write('''ERROR: indelible not found. Needs to be installed to simulate insertations
                                        and deletions. IGNORING indel parameters and continuaing simulation\n''')
                    self.mut_genomes_no_indels()
                else:
                    self.mut_genomes_indels()
            else:
                self.mut_genomes_no_indels()
        if coverage is None:
            coverarg = self.get_arg('cov')
        else:
            coverarg = coverage
        cov = {}
        if os.path.isfile(coverarg):
            with open(coverarg) as infile:
                for lin in infile:
                    seqnam = lin.split(',')[0]
                    try:
                        assert seqnam in self.seqnames
                    except:
                        sys.stderr.write("name {} in coverage file not found in tree\n".format(seqnam))
                        self._exit_handler()
                    cov[seqnam] = int(lin.split(',')[1])
            try:
                assert set(cov.keys()) == set(self.seqnames)
            except:
                sys.stderr.write("some tips missing from coverage file: {}\n".format(set(self.seqnames) - cov.keys()))
                self._exit_handler()
        else:
            for seqnam in self.seqnames:
                cov[seqnam] = int(coverarg)
        if not os.path.isdir("{}/fastq".format(self.outd)):
            os.mkdir("{}/fastq".format(self.outd))
        sys.stdout.write("coverage is {}\n".format(self.config['coverage']))
        if 'read_length' in self.config:
            read_length = self.config['read_length']
        else:
            read_length = 150
        sys.stdout.write("read length is {}\n".format(read_length))
        if 'fragment_size' in self.config:
            fragment_size = self.config['fragment_size']
        else:
            fragment_size = 350
        sys.stdout.write("fragment size is {}\n".format(fragment_size))
        if 'stdev_frag_size' in self.config:
            stdev_frag_size = self.config['stdev_frag_size']
        else:
            stdev_frag_size = 130
        sys.stdout.write("stdev of frag size is {}\n".format(stdev_frag_size))
        for seq in self.seqnames:
            sys.stdout.write("Generating reads for {}\n".format(seq))
            if not os.path.isdir("{}/fastq/{}{}".format(self.outd, self.prefix, seq)):
                os.mkdir("{}/fastq/{}{}".format(self.outd, self.prefix, seq))
            if  self.get_arg('error_model1') and self.get_arg('error_model2'):
                assert os.path.exists(self.get_arg('error_model1'))
                assert os.path.exists(self.get_arg('error_model2'))
                artparam = ['art_illumina',
                            '-1', self.get_arg('error_model1'),
                            '-2', self.get_arg('error_model2'),
                            '-na', #Don't output alignment file
                            '-p', #for paired end reads
                            '-i', '{}/fasta_files/{}{}.fasta'.format(self.outd, self.prefix, seq),
                            '-l', '{}'.format(read_length),
                            '-f', str(cov[seq]),
                            '-m', '{}'.format(fragment_size),
                            '-s', '{}'.format(stdev_frag_size),
                            '-o', '{}/fastq/{}{}/{}{}_'.format(self.outd,
                                                               self.prefix,
                                                               seq,
                                                               self.prefix,
                                                               seq)]
            else:
                artparam = ['art_illumina',
                            '-p', #for paired end reads
                            '-na', #Don't output alignment file
                            '-i', '{}/fasta_files/{}{}.fasta'.format(self.outd, self.prefix, seq),
                            '-l', '{}'.format(read_length),
                            '-f', str(cov[seq]),
                            '-m', '{}'.format(fragment_size),
                            '-s', '{}'.format(stdev_frag_size),
                            '-o', '{}/fastq/{}{}/{}{}_'.format(self.outd,
                                                               self.prefix,
                                                               seq,
                                                               self.prefix,
                                                               seq)]
            call(artparam, stdout=open('{}/art_log'.format(self.outd), 'w'), stderr=open('{}/art_log'.format(self.outd), 'a'))
          #  print("called {}".format(" ".join(artparam)))
            assert os.path.exists('{}/fastq/{}{}/{}{}_1.fq'.format(self.outd, self.prefix, seq, self.prefix, seq))
            gzippar = ['gzip',
                       '-f',
                       '{}/fastq/{}{}/{}{}_1.fq'.format(self.outd,
                                                        self.prefix,
                                                        seq,
                                                        self.prefix,
                                                        seq),
                       '{}/fastq/{}{}/{}{}_2.fq'.format(self.outd,
                                                        self.prefix,
                                                        seq,
                                                        self.prefix,
                                                        seq)]
            call(gzippar)
        sys.stdout.write("TreeToReads completed successfully!\n")




def write_indelible_controlfile(outputdir, ratemat, freqmat, indelmodel, indelrate, tree, seqlen, seed):
    """Writes a control file for indelible to run
         indelible GTR is specified as ct', 'at', 'gt', 'ac', 'cg', 'ag' =1
         and indelible base frequencies are specified as 
         pi_T, pi_C, pi_A, pi_G"""
    reorder_freqs = ' '.join([freqmat[base] for base in ['T','C','A','G']])
    rescale_rates = ' '.join([str(float(ratemat[rate])/float(ratemat['ag'])) for rate in  ['ct', 'at', 'gt', 'ac', 'cg', 'ag']])
    fi = open("{}/control.txt".format(outputdir), 'w')
    fi.write("[TYPE] NUCLEOTIDE 1 \n")
    fi.write("[SETTINGS]\n")
    fi.write("[output] FASTA \n")
    fi.write("[randomseed]  {}\n".format(seed))
    fi.write("[MODEL] TTRm \n")
    fi.write("[submodel] GTR {} \n".format(rescale_rates))
    fi.write("[indelmodel] {} \n".format(indelmodel))
    fi.write("[indelrate]  {} \n".format(indelrate))
    fi.write("[statefreq]  {} \n".format(reorder_freqs))
    fi.write("[TREE] TTR {};\n".format(tree))
    fi.write("[PARTITIONS] partitionname\n")
    fi.write("[TTR TTRm {}]\n".format(int(seqlen*1.2)))
    fi.write("[EVOLVE] partitionname 1 TTRindelible\n")
    fi.close()


def run_indelible(outputdir):
    """Calls indelible, and returns to working directory"""
    cwd = os.getcwd()
    os.chdir(outputdir)
    call(['indelible', "control.txt", ">", "indelible.log"])
    os.chdir(cwd)


def read_indelible_aln(ttrobj):
    """Pull steh locations of insertaions and deletions from indelible output.
    Insertion locs are in terms of the original sequence length,
    but deletions are in terms of alignment length"""
    insertionlocs = {}
    insertionlocs_aln = set()
    insertions = {}
    deletions = {}
    base_name = ttrobj.get_arg('base_name')
    indel_aln = open("{}/TTRindelible_TRUE.fas".format(ttrobj.outd))
    for lin in indel_aln:
        if lin.startswith(">"):
            seqname = lin.strip(">").strip().strip("'")
            assert seqname in ttrobj.seqnames
        elif seqname == base_name:
            ref_genome_i = 0
            for i, char in enumerate(lin):
                if char == '-':
                    insertionlocs_aln.add(i)
                    if not insertionlocs.get(ref_genome_i):
                        insertionlocs[ref_genome_i] = set()
                        insertionlocs[ref_genome_i].add(i)
                    else:
                        insertionlocs[ref_genome_i].add(i)
                else:
                    ref_genome_i += 1
                if ref_genome_i == ttrobj.genlen:
                    alignment_length = i
                    sys.stdout.write("Base genome length is  {} and alignment length will be {}\n".format(ttrobj.genlen, i))
                    break
    indel_aln = open("{}/TTRindelible_TRUE.fas".format(ttrobj.outd))
    del_locs = {}
    for lin in indel_aln:
        if lin.startswith(">"):
            seqname = lin.strip(">").strip().strip("'")
            del_locs[seqname]=set()
            #assert seqname in ttrobj.seqnames
        elif seqname and seqname != base_name:
            insertions[seqname] = {}
            deletions[seqname] = set()
            for i, char in enumerate(lin):
                if i in insertionlocs_aln:
                    insertions[seqname][i] = char
                elif char == '-':
                    deletions[seqname].add(i) ##
                    if seqname in del_locs:
                        del_locs[seqname].add(i)
                    else:
                        del_locs[seqname] = set()
                        del_locs[seqname].add(i)
                if i >= alignment_length:
                    break
            seqname = None
        deletions[base_name] = {}
  #  print del_locs
    total_dellocs = set()
    for seqname in del_locs:
        sub_del_locs = list(del_locs[seqname])
        sub_del_locs.sort()
        tmpout.write("\n"+seqname+"ln893\n")
        tmpout.write(", ".join([str(i) for i in sub_del_locs]))
        tmpout.write("\n")
   #     print sub_del_locs
        sub_deletionlocs = get_sub_list(sub_del_locs)
    #    print "sublist looks like"
     #   print sub_deletionlocs
        tmpout.write(", ".join([str(i) for i in sub_deletionlocs]))
        fixed_dels = [tuple(i) for i in sub_deletionlocs if len(i) >= 1]
        deletions[seqname] = fixed_dels
        for deletion in fixed_dels:
            total_dellocs.add(deletion)
    list_of_deletions = list(total_dellocs)
    list_of_deletions = sorted(list_of_deletions,key=itemgetter(0))
   # print "list of deeltion"
   # print list_of_deletions
    doublehit = set()
    del_DH = []
    translate_deletions = {}
    deletionlocs = []
    for i, dele in enumerate(list_of_deletions):
        for base in dele:
            if base in deletionlocs:
               print "DOUBLEHIT deletion {}".format(base)
               del_DH.append(base)
            else:
                deletionlocs.append(base)
    return insertions, deletions, insertionlocs, list_of_deletions

def split_list(n):
    """will return the list index for sequential deletions"""
    return [(x+1) for x, y in zip(n, n[1:]) if y-x != 1]

def get_sub_list(my_list):
    """will split the list based on the index for splits"""
    my_index = split_list(my_list)
    output = list()
    prev = 0
    for index in my_index:
        new_list = [x for x in my_list[prev:] if x < index]
        output.append(new_list)
        prev += len(new_list)
    output.append([x for x in my_list[prev:]])
    return output

def write_vcf(ttrobj):
    """Writes out simulated mutations and indels as a vcf file, with anchor genome as reference."""
    fi = open("{}/sim.vcf".format(ttrobj.outd), 'w')
    fi.write("##fileformat=VCFv4.0\n")
    fi.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fi.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format('\t'.join(ttrobj.seqnames)))
    mutlocs = ttrobj.vcf_dict.keys()
    mutlocs.sort()
#    assert mutlocs == ttrobj.mutlocs
   # contig = 0
    for loc in mutlocs:
        #adjusted_loc = loc - ttrobj.contig_breaks[contig]
        #if ttrobj.contig_breaks[contig+1] < loc:
        #    contig += 1
        #adjusted_loc = loc - ttrobj.contig_breaks[contig]
        #contig_name = ttrobj.contig_names[contig]
        #ttrobj.mut_trans[loc] = (contig_name, str(adjusted_loc))
        contig_name, adjusted_loc = ttrobj.mut_trans[loc]
        assert set(ttrobj.vcf_dict[loc].keys()) == set(ttrobj.seqnames)
        refbase = ttrobj.vcf_dict[loc][ttrobj.base_name]
        base_calls = [ttrobj.vcf_dict[loc][seq] for seq in ttrobj.seqnames]
        for i, nuc in enumerate(base_calls):
            if '-' in nuc:
                base_calls[i] = nuc.replace('-', '.')
            if nuc.replace('-', '') == refbase:
                base_calls[i] = refbase
        altbase = set(base_calls) - set([refbase])
        trans = {refbase:'0'}
        for i, base in enumerate(altbase):
            trans[base] = str(i+1)
        variants = [trans[base] for base in base_calls]
        fi.write('''{chrm}\t{loc}\t.\t{refbase}\t{altbase}\t40\tPASS\t.\tGT\t{vars}\n'''.format(chrm=contig_name,
                                                       loc=int(adjusted_loc)+1,
                                                       refbase=refbase,
                                                       altbase=",".join(altbase),
                                                       vars='\t'.join(variants)))
    fi.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Tree to Reads - A python script to to read a tree,
                    resolve polytomes, generate mutations and simulate reads.""",
        epilog="""Still in development - email ejmctavish@gmail.com with questions, suggestions, issues etc.""")
    parser.add_argument("config_file", type=str, help="configuration file path. Required, defaults to seqsim.cfg")
    parser.add_argument('-V', '--version',
                        action='version',
                        version='Tree to reads version {}'.format(VERSION))
    args = parser.parse_args()
    ttr = TreeToReads(configfi=args.config_file, main=1)


