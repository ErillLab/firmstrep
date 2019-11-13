# -*- coding: utf-8 -*-
"""
Reads in a comma-separated file containing the genbank accession and other
details about phage genomes. Downloads the corresponding genome files in
both GenBank (full) and FASTA format, stores them in the corresponding 
'genbank' and 'fasta' directories under 'genomes'.
"""

########################### IMPORTS #######################################
#import json to write out json file
import json

#system library to time out
import time

#import system lib to check dir contents
from os import listdir

#CSV library to read off the csv file
import csv

#imports the Entrez and SeqIO modules from the Bio (biopython) library
#Bio is the official name for biopython
#we also set the basic entrez parameters (e.g. email address)
from Bio import Entrez
Entrez.email ="ivan.erill@gmail.com"

#seqIO for GBK to FASTA conversion
from Bio import SeqIO

#bioutils to compute %GC
from Bio.SeqUtils import GC

#import genome parsing from BINFUtils
import Genome_parsing as GP

#genomes path
gnome_path = '../genomes'

#data path
data_path = '../src_data'

#out data path
out_path = '../derived_data'

phages={}

#read in the csv file into the phage dictionary
with open(data_path+'/table.csv', 'rb') as csvfile:
    csvobj = csv.DictReader(csvfile)
    for row in csvobj:
        phages[row['Name']]=row
        
        
#for each genome, download GBK, and save using phage name
for phage, attributes in phages.iteritems():
    prophage=False
    dir_contents=listdir(gnome_path+'/genbank/')
    #for each phage download accession to file name.gb, name.fas
    print "Processing: " + attributes['Accession']
    #check whether the genome is already in
    if not((attributes['Accession']+'.gb') in dir_contents):
        print "--> Downloading: " + attributes['Accession']
        #retrieve genome record
        net_handle = Entrez.efetch(db="nuccore",id=attributes['Accession'],
                                   rettype='gbwithparts', retmode="txt")
        gnome_record=net_handle.read()
        #write the record to file
        out_handle = open(gnome_path+'/genbank/'+attributes['Accession']+".gb", "w")
        out_handle.write(gnome_record)
        time.sleep(1.5)
        out_handle.close()
        
        out_handle.close()
        net_handle.close()

    #generate FASTA sequence file from genbank one
    SeqIO.convert(gnome_path+'/genbank/'+attributes['Accession']+".gb",
                  "genbank", 
                  gnome_path+'/fasta/'+attributes['Accession']+".fas", 
                  "fasta")        
    
    #read up genbank file and update fields        
    gbkrecord = SeqIO.read(gnome_path+'/genbank/'+attributes['Accession']
                           +'.gb', 'genbank')
    
    #generate AA FASTA file from GenBank record
    GP.GBK2mFASTA(gbkrecord, gnome_path+'/fasta_aa/'\
                  +attributes['Accession']+".fas", f_type='CDS', s_type='AA')
    
    #deal with prophages
    prophage=False
    if len(gbkrecord.features[0])>1500000:
        prophage=True
        attributes['Prophage']=1
    else:
        attributes['Prophage']=0

    #get host
    attributes['Plasmid']=0
    if 'host' in gbkrecord.features[0].qualifiers:
        attributes['Host']=gbkrecord.features[0].qualifiers['host'][0]
    else:
        if 'lab_host' in gbkrecord.features[0].qualifiers:
            attributes['Host']=gbkrecord.features[0].qualifiers['lab_host'][0]
        else:
            #deal with plasmid
            if 'plasmid' in gbkrecord.features[0].qualifiers:
                attributes['Host']=gbkrecord.features[0].qualifiers['organism'][0]
                attributes['Plasmid']=1
            else:
                attributes['Host']=None
       
    #get morphotype
    if 'taxonomy' in gbkrecord.annotations:
        if attributes['Plasmid']==0:
            #if detailed taxonomy is provided
            if 'dsDNA viruses, no RNA stage' in gbkrecord.annotations['taxonomy']:
                #get position for dsDNA virus term
                cnt=0
                while gbkrecord.annotations['taxonomy'][cnt]!='dsDNA viruses, no RNA stage':
                    cnt=cnt+1
        
                attributes['Order']=gbkrecord.annotations['taxonomy'][cnt+1]
                if (cnt+2)<len(gbkrecord.annotations['taxonomy']):
                    attributes['Family']=gbkrecord.annotations['taxonomy'][cnt+2]
                else:
                   attributes['Family']=None 
            else:
                print "No taxonomy reported for: ", phage
                attributes['Order']=None
                attributes['Family']=None
        else:
            attributes['Order']=gbkrecord.annotations['taxonomy'][-4]
            attributes['Family']=gbkrecord.annotations['taxonomy'][-3]
    else:
        print "No detailed taxonomy reported for: ", phage
        attributes['Order']=None
        attributes['Family']=None
            
    #get country
    if 'country' in gbkrecord.features[0].qualifiers:
        attributes['Country']=gbkrecord.features[0].qualifiers['country'][0].split(':')
    else:
        attributes['Country']=None

    #get coordinates
    if 'lat_lon' in gbkrecord.features[0].qualifiers:
        attributes['Coordinates']=gbkrecord.features[0].qualifiers['lat_lon'][0].split()
    else:
        attributes['Coordinates']=[None,None,None,None]
        
         
    #get genome size
    attributes['Genome_length']=len(gbkrecord.features[0])
    
    #get %GC content
    attributes['GC_content']=GC(gbkrecord.seq)
    
#    print 'Order: ', attributes['Order']
#    print 'Family: ', attributes['Family']
#    print 'Host: ', attributes['Host']
#    print 'eHost: ', attributes['eHost']
#    print 'Country: ', attributes['Country']
#    print 'Coordinates: ', attributes['Coordinates']
#    print 'GroupA: ', attributes['GroupA']
#    print 'GroupB: ', attributes['GroupB']
#    print 'Genome_length: ', attributes['Genome_length']
#    print 'GC_content: ', attributes['GC_content']
#    print 'Cluster :', attributes['Cluster']
#    print 'Subcluster :', attributes['Subcluster']
#    print 'Plasmid :', attributes['Plasmid']


with open(out_path+'/Processed_table.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    #write dict keys to first csv line
    writer.writerow(phages.itervalues().next().keys())
    #write data lines
    for key, value in phages.iteritems():
        writer.writerow(value.values())    

phage_list=[]
for key, value in phages.iteritems():
    phage_list.append(value)
    
with open(out_path+'/Phage_data.json', 'wb') as jsonfile:
    json.dump(phage_list, jsonfile, sort_keys=True, indent=4)
        

    
#        #update the json accession list to contain RefSeq ID
#        new_access.append(accession)
#        out_handle.close()
#        net_handle.close()
#    #if the genome was already in, do not update json with empty accessions
#    else:
#        updatejson=False
#
##update json only if we did download genome and obtained RefSeq accessions
#if updatejson:
#    #update the json dict
#    orgn['RefSeq_accessions']=new_access
#    #add bool to trace that record has been mapped to RefSeq in the json file
#    orgn['RefSeq_mapped']=True
#    
#    #save json dict list
#    with open(data_path+'/'+'dataset.json', 'w') as outfile:  
#        json.dump(dataset, outfile, indent=4, sort_keys=True, separators=(',', ':'))    
#        
##    
#    for k,v in csvobj:
#        phages[k]=v
#        
#        
#for row in reader:
#   k, v = row
#   d[k] = v
