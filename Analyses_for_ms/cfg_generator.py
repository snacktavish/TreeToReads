#Generate (and run?) the full set of CFG files to create data for anayses


#Nameformat Coverage, number, treefile, run e.g. C1_N200_T1_R0.cfg



'''Avg. coverage 1x
Avg. coverage 10x
Avg coverage 50X
Seq. Div. Reality (200) (decide on exp param)
Seq. Div. 1% (40,000) (decide on exp param)
Seq. Div. 5%  (200,000) (decide on exp param)
OB node: 1 SNP on branch leading to outbreak (treefile1)
OB node: 5x SNP (treefile2)
OB node: 10x SNP (treefile3)
Ref. Reality (internal to tree)
Ref. 1%
Ref. 5%
'''

argdict = {'treepath':'treefile_path', 'nsnp':'number_of_snps', 'anchor_name':'anchor_name', 'genome':'anchor_genome_path', 'ratmat':'rate_matrix', 'freqmat':'freq_matrix', 'shape':'shape', 'errmod1':'error_model1', 'errmod2':'error_model2', 'cov':'coverage', 'outd':'output_dir'}



default = {'treepath':'tree_AlleleCounts.ML.tre', 'nsnp': 20 , 'anchor_name':'CFSAN000189', 'genome':'../example/barref.fasta', 'ratmat':'0.441,5.435,0.717,0.020,6.001,1.000',\
 'freqmat':'0.19,0.31,0.29,0.22', 'shape':'50', 'errmod1':'../example/ErrprofR1.txt', 'errmod2':'../example/ErrprofR2.txt', 'cov':20, 'outd':'output_dir'}


def makecfg(params):
    fi=open('{}.cfg'.format(params['id']),'w')
    for item in argdict:
    	fi.write("{} = {}\n".format(argdict[item],params[item]))
    fi.close()


run=0
for coverage in [1,10,50]:
    for number in [200,40000,100000]:
        for tree in ['T1','T2','T3']:
            params= dict(default)
            params['id']='C{}_N{}_{}_R{}'.format(coverage,number,tree,run)
            params['outd']=params['id']
            params['cov']=coverage
            params['nsnp']=number
            params['treepath']='{}.tre'.format(tree)
            makecfg(params)




