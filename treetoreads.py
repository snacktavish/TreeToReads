#!/usr/bin/env python
"""Tree to Reads - A python script to to read a tree,
    resolve polytomies, generate mutations and simulate reads."""
import random
import os
import sys
import argparse
from subprocess import call
import dendropy


VERSION = "0.0.3"

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
            if (self.no_art == 0) and (int(self.config['coverage']) != 0):
                self.run_art()
            else:
                sys.stdout.write("simulating genomes but not reads\n")
                self.mut_genomes()
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
        if call(['which', 'samtools'], stdout=open('/dev/null', 'w')) == 1:
            sys.stderr.write('''ERROR: samtools needs to be installed to output sorted bam files.
                             samtools not found
                             TTR will generate SAM files (and may take up a lot of disk space!) \n''')
            self.no_sam = 1
        else:
            self.no_sam = 0

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
                     'error_model2']
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
                        'ratmat':'rate_matrix',
                        'freqmat':'freq_matrix',
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
        if not len(self.get_arg('ratmat').split(',')) == 6:
            sys.stderr.write('''{} values in rate matrix, there should be 6.
                              Exiting.\n'''.format(len(self.get_arg('ratmat').split(','))))
            self._exit_handler()
        if not len(self.get_arg('freqmat').split(',')) == 4:
            sys.stderr.write('''{} values in freq matrix, there should be 4.
                              Exiting.\n'''.format(len(self.get_arg('freqmat').split(','))))
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
        #import tree from path
        if dendropy.__version__.startswith('4'):
            taxa = dendropy.TaxonNamespace()
            try:
                tree = dendropy.Tree.get_from_path(self.get_arg('treepath'), self.treetype, taxon_namespace=taxa)
            except:
                sys.stderr.write("Problems reading the tree - is it in proper newick or nexus format?\n")
        else:
            taxa = dendropy.TaxonSet()
            try:
                tree = dendropy.Tree.get_from_path(self.get_arg('treepath'), self.treetype, taxon_set=taxa)
            except:
                sys.stderr.write("Problems reading the tree - is it in proper newick or nexus format?\n")
        if tree.length() == 0:
            sys.stderr.write("TTR requires branch lengths. Branch lengths appear to be missing (treelength = 0). Exiting.\n")
            self._exit_handler()
        self.seqnames = taxa.labels()
        if not self.get_arg('base_name') in self.seqnames:
            sys.stderr.write("base genome name {} is not in tree. Exiting.\n".format(self.get_arg('base_name')))
            self._exit_handler()
        tree.resolve_polytomies()
        self.seqgen_scaler = float((1.0/tree.length()) * 10)
#        sys.stdout.write("scaler is {}\n".format(self.seqgen_scaler))
        if tree.length() >= 1:
            sys.stderr.write("WARNING: Tree length is high- scale down tree or expect high multiple hits/homoplasy!\n")
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
        with open(self.get_arg('genome'), 'r') as in_file:
            for line in in_file:
                line = line.strip()
                if line.startswith('>'):
                    contigs += 1
                    self.contig_breaks.append(self.genlen)
                else:
                    self.genlen += len(line)
                    if not set(line.upper()).issubset(set(['A', 'T', 'G', 'C','N'])):
                        sys.stderr.write('''Your genome appears to have characters other than ATGC,
                                         such as: {} Please check your input genome.\n'''.format(set(line)))
        if self.nsnp > self.genlen:
            sys.stderr.write('''number of variables sites {}
                             is higher than the length
                             of the contig or geonme {}.
                             Exiting\n'''.format(self.nsnp, self.genlen))
            self._exit_handler()
        sys.stdout.write("{} contigs\n".format(len(self.contig_breaks)))
        sys.stdout.write("Genome has {} bases\n".format(self.genlen))

    def generate_varsites(self):
        """Runs seqgen to generate variable sites on tree"""
        self._vargen = 1
        if not self._treeread:
            self.read_tree()
        self.simloc = "{}/seqs_sim.txt".format(self.outd)
        lenseqgen = self.seqgen_scaler*self.nsnp
        if self.shape:
            seqgenpar = ['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR',
                         '-a{}'.format(self.shape), '-r{}'.format(self.get_arg('ratmat')),
                         '-f{}'.format(self.get_arg('freqmat')), '-or']
        else:
            seqgenpar = ['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR',
                         '-r{}'.format(self.get_arg('ratmat')),
                         '-f{}'.format(self.get_arg('freqmat')), '-or']
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
        if add == False:
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
            sys.stderr.write(self.seqnames)
            sys.stderr.write(simseqs.keys())
            self._exit_handler()
        ref = simseqs[self.get_arg('base_name')]
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
        sys.stderr.write('simulating additional sites- consider scaling up input tree branch lengths\n')
        if self._siteread != 1:
            self.read_sites()
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
            ran = random.sample(range(self.genlen), 1+(self.nsnp-len(rands)))
            rands = rands | set(ran)
        for site in rands:
            fi.write(str(site)+'\n')
        self.mutlocs = rands
        self._mutlocs = 1

    def assign_sites(self):
        """Pulls columns for each mutation"""
        if not self._mutlocs:
            self.select_mutsites()
        if not self._siteread:
            self.read_varsites()
        self.mut_genos = {}
        patnuc = {}
        ri = 0
        patnuc['A'] = 0
        patnuc['G'] = 0
        patnuc['T'] = 0
        patnuc['C'] = 0
        self.snpdic = {}
        #This next section is to index mutations, and re-simulate if there are not enough of a base
        with open(self.get_arg('genome'), 'r') as in_file:
            for line in in_file:
                line = line.strip()
                if line.startswith('>'):
                    pass
                else:
                    for nuc in line:
                        ri += 1
                        if ri in self.mutlocs:
                            patnuc[nuc] += 1
                            self.snpdic[ri] = patnuc[nuc]
                            if len(self.sitepatts[nuc]) < patnuc[nuc]:
                                self.add_varsites()
        sys.stderr.write('''{} variable sites with more than 2 bases out of {} sites.
        Scale down tree length if this is too high.\n'''.format(self.trip_hit, self.var_site))


    def mut_genomes(self):
        """Writes out the simulated genomes with mutations"""
        self.assign_sites()
        matout = open("{}/var_site_matrix".format(self.outd), 'w')
        for seq in self.seqnames:
            self.mut_genos[seq] = []
            sys.stdout.write("writing genome for {}\n".format(seq))
            if not os.path.isdir("{}/fasta_files".format(self.outd)):
                os.mkdir("{}/fasta_files".format(self.outd))
            genout = open("{}/fasta_files/{}{}.fasta".format(self.outd, self.prefix, seq), 'w')
            ii = 0
            with open(self.get_arg('genome'), 'r') as in_file:
                for line in in_file:
                    if line.startswith('>'):
                        if ii > 0:
                            genout.write('\n')
                        genout.write(line.strip()+"_"+self.prefix+seq)
                        ii = 0
                    else:
                        line = line.strip()
                        for nuc in line:
                            if ii%70 == 0:
                                genout.write('\n')
                            ii += 1
                            if ii in self.mutlocs:
                                if nuc == 'N':
                                    genout.write('N')
                                else:
                                    patt = self.sitepatts[nuc][self.snpdic[ii]]
                                    genout.write(patt[seq])
                                    self.mut_genos[seq].append(patt[seq])
                                    matout.write("{} {} {}\n".format(seq, patt[seq], ii))
                            else:
                                genout.write(nuc)
            genout.write('\n')
            genout.write('\n')
            genout.close()
        matout.close()
        self._genmut = 1
        sys.stdout.write("Mutated genomes\n")

    def run_art(self):
        """Runs ART to simulate reads from the simulated genomes"""
        if not self._genmut:
            self.mut_genomes()
        if not os.path.isdir("{}/fastq".format(self.outd)):
            os.mkdir("{}/fastq".format(self.outd))
        if not self._genmut:
            self.mut_genomes()
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
                            '-sam', #output a sam file
                            '-p', #for paired end reads
                            '-i', '{}/fasta_files/{}{}.fasta'.format(self.outd, self.prefix, seq),
                            '-l', '{}'.format(read_length),
                            '-f', self.get_arg('cov'),
                            '-m', '{}'.format(fragment_size),
                            '-s', '{}'.format(stdev_frag_size),
                            '-o', '{}/fastq/{}{}/{}{}_'.format(self.get_arg('outd'),
                                                               self.prefix,
                                                               seq,
                                                               self.prefix,
                                                               seq)]
            else:
                artparam = ['art_illumina',
                            '-p', #for paired end reads
                            '-na', #Don't output alignment file
                            '-sam', #output a sam file
                            '-i', '{}/fasta_files/{}{}.fasta'.format(self.outd, self.prefix, seq),
                            '-l', '{}'.format(read_length),
                            '-f', self.get_arg('cov'),
                            '-m', '{}'.format(fragment_size),
                            '-s', '{}'.format(stdev_frag_size),
                            '-o', '{}/fastq/{}{}/{}{}_'.format(self.get_arg('outd'),
                                                               self.prefix,
                                                               seq,
                                                               self.prefix,
                                                               seq)]
            call(artparam, stdout=open('{}/art_log'.format(self.outd), 'w'))
            assert os.path.exists('{}/fastq/{}{}/{}{}_1.fq'.format(self.get_arg('outd'), self.prefix, seq, self.prefix, seq))
            gzippar = ['gzip',
                       '-f',
                       '{}/fastq/{}{}/{}{}_1.fq'.format(self.get_arg('outd'),
                                                        self.prefix,
                                                        seq,
                                                        self.prefix,
                                                        seq),
                       '{}/fastq/{}{}/{}{}_2.fq'.format(self.get_arg('outd'),
                                                        self.prefix,
                                                        seq,
                                                        self.prefix,
                                                        seq)]
            call(gzippar)
            if not self.no_sam:
                # LK - sam => bam
                sam = '{}/fastq/{}{}/{}{}_.sam'.format(self.get_arg('outd'),
                                                       self.prefix,
                                                       seq,
                                                       self.prefix,
                                                       seq)
                bam = '{}/fastq/{}{}/{}{}.bam'.format(self.get_arg('outd'),
                                                      self.prefix,
                                                      seq,
                                                      self.prefix,
                                                      seq)
                sortedBamPrefix = '{}/fastq/{}{}/{}{}.sorted'.format(self.get_arg('outd'),
                                                                     self.prefix,
                                                                     seq,
                                                                     self.prefix,
                                                                     seq)
                sortedBam = '{}/fastq/{}{}/{}{}.sorted.bam'.format(self.get_arg('outd'),
                                                                   self.prefix,
                                                                   seq,
                                                                   self.prefix,
                                                                   seq)
                # convert sam to bam
                call(['samtools', 'view', '-bS', '-o', bam, sam])
                # sort the bam
                call(['samtools', 'sort', bam, sortedBamPrefix])
                # index the bam
                call(['samtools', 'index', sortedBam])
                # remove intermediate files
                os.remove(sam)
                os.remove(bam)
        sys.stdout.write("TreeToReads completed successfully!\n")




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
