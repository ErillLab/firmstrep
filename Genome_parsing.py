# -*- coding: utf-8 -*-

"""
Created on Sept 17 2018

Functions to extract specific records from genome files and save them in 
specific formats.

@author: ivanerill
"""
#imports the Entrez module from the Bio (biopython) package
#Bio is the official name for biopython
#we also set the basic entrez parameters (e.g. email address)
from Bio import Entrez, SeqIO
Entrez.email ="ivan.erill@gmail.com"

def NCBIreadGBK(accession):
    """Gets an accession or GI number and downloads the GenBank record.
       Parses it with SeqIO and returns the sequence object.
    """
    net_handle = Entrez.efetch(db="nuccore",id=str(accession),
                                   rettype='gbwithparts', retmode="txt")
    gnome_record=SeqIO.read(net_handle, "genbank")
    net_handle.close()
    return gnome_record

def readGBK(filename):
    """Reads a GBK file and returns the object handle"""
    gnome_record = SeqIO.read(filename, "genbank")
    return gnome_record

def GBK2mFASTA(genomeobject, outfilename, f_type='CDS', s_type='AA'):
    """Processes a GBK record object, looking for genome features and saves 
       their nucleotide or amino acid sequence in multiline FASTA format.
       
       Inputs:
       - genomeobject
           a parsed GBK genome object, as done by SeqIO.read (e.g. readGBK)
       - outfilename
           the output file name. script does NOT add '.fas' extension
       - f_type
           a feature type
           currently defined: 'CDS', 'tRNA', 'rRNA'
           feel free to add more in 'feat_types' below           
       - s_type
           a sequence type
           currently only DNA and AA supported
    """
    
    #define acceptable feature types
    feat_types = ['CDS', 'tRNA', 'rRNA']
    
    
    #check feature type
    if f_type not in feat_types:
        print 'WARNING: Unrecognized feature type <' + f_type \
                  + '>. Execution stopped.'
        return(-1)
                  
    #check that AA is used only for CDS; otherwise switch to DNA
    if s_type == 'AA':
        if f_type != 'CDS':
            s_type = 'DNA'
            print 'WARNING: Sequence type for feature type <' + f_type \
                  + '> has been switched to <DNA>'

              
    with open(outfilename,"w") as out_handle:
        ft_cnt=1
        #iterate features
        for feat in genomeobject.features:
            #get type features
            if feat.type == f_type:
                out_handle.write('>')
                out_handle.write(genomeobject.id + ' | ')
                if f_type == 'CDS':
                    if 'protein_id' in feat.qualifiers:
                        out_handle.write(feat.qualifiers['protein_id'][0] + ' | ')
                if 'product' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['product'][0] + ' | ')
                if 'locus_tag' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['locus_tag'][0] + ' | ')
                else:
                    out_handle.write('NO_LOCUS_TAG_' + str(ft_cnt) + ' | ')
                out_handle.write(str(feat.location.nofuzzy_start) + '-')
                out_handle.write(str(feat.location.nofuzzy_end) + ' : ')
                out_handle.write(str(feat.location.strand))
                if 'old_locus_tag' in feat.qualifiers:
                    out_handle.write('| ['+feat.qualifiers['old_locus_tag'][0]+']')
                out_handle.write('\n')
                if s_type == 'AA':
                    out_handle.write(str(feat.qualifiers['translation'][0]))
                else:
                    out_handle.write(str(feat.location.extract(genomeobject).seq))
                out_handle.write('\n')   
                ft_cnt=ft_cnt+1
            
    return(0)