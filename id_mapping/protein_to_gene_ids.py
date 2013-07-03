# -*- coding: utf-8 -*-

def my_efetch(db, **keywds):
    """Fetches Entrez results which are returned as a handle.

    EFetch retrieves records in the requested format from a list of one or
    more UIs or from user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

    Return a handle to the results.

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.efetch(db="nucleotide", id="57240072", rettype="gb", retmode="text")
    >>> print handle.readline().strip()
    LOCUS       AY851612                 892 bp    DNA     linear   PLN 10-APR-2007
    >>> handle.close()

    Warning: The NCBI changed the default retmode in Feb 2012, so many
    databases which previously returned text output now give XML.
    """
    cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    variables = {'db': db}
    keywords = keywds
    post=False
    if "post" in keywds: post=keywds["post"]
    if "id" in keywds and isinstance(keywds["id"], list):
        #Fix for NCBI change (probably part of EFetch 2,0, Feb 2012) where
        #a list of ID strings now gives HTTP Error 500: Internal server error
        #This was turned into ...&id=22307645&id=22303114&... which used to work
        #while now the NCBI appear to insist on ...&id=22301129,22299544,...
        keywords = keywds.copy()  # Don't alter input dict!
        keywords["id"] = ",".join(keywds["id"])
    variables.update(keywords)
    return Entrez._open(cgi, variables, post)


from Bio import Entrez
import re
query_db='protein'
Entrez.efetch=my_efetch


def acc_to_gi(id_list, rettype):
    results=[]
    if(len(id_list)):
        results_handle=Entrez.efetch(db=query_db, rettype=rettype, retmode="text", post=True, id=id_list)
        result_list=results_handle.readlines()
        try:
            assert len(id_list)==len(result_list)
        except(AssertionError):
            print "couldn't convert all accessions to gi's. there may be a problem"
        for idx, val in enumerate(result_list):
            if rettype=="acc" and '.' in val:
                result_list[idx]=val.replace(val[val.find('.'):],'')
            else:
                result_list[idx]=val.strip()
        return result_list


## id list should be accessions
def get_gene_ids(id_list):
    assert re.search(r'[a-zA-Z]{2}_[0-9]+',id_list[0])
    results=id_list[:]
    for idx, val in enumerate(results): results[idx]=(results[idx],'-')
    order_ids={}
    for idx, val in enumerate(id_list): order_ids[val]=idx
    if(len(id_list)):
        results_handle=Entrez.efetch(db=query_db, rettype="gp", retmode="xml", post=True, id=id_list)
        result_list=Entrez.read(results_handle)
        try:
            assert len(id_list)==len(result_list)
        except(AssertionError):
            print "records retrieved: "+str(len(result_list))+" not equal to IDs: "+str(len(id_list))
        for rec in result_list:
            cur_acc=rec['GBSeq_primary-accession']
            if cur_acc in order_ids:
                idx=order_ids[cur_acc]
                if 'GBSeq_feature-table' in rec:
                    for feat in rec['GBSeq_feature-table']:
                        if 'GBFeature_quals' in feat:
                            for qual in feat['GBFeature_quals']:
                                if 'GBQualifier_name' in qual and qual['GBQualifier_name']=='db_xref' \
                                and 'GBQualifier_value' in qual and qual['GBQualifier_value'].startswith('GeneID:'):
                                    value= qual['GBQualifier_value'].split(':')[-1]
                                    results[idx]= (results[idx][0], value)
        return results


def get_updated(gi_list, rettype='acc', remove_old=False):
    results=gi_list[:]
    order_ids={}
    for idx, val in enumerate(gi_list): order_ids[val]=idx
    if rettype=='acc':
        rfield="Caption"
    else: rfield="Gi"
    post_handle=Entrez.epost(db='protein', id=','.join(gi_list))
    post_result=Entrez.read(post_handle)
    summary_handle=Entrez.esummary(db=query_db, query_key=post_result['QueryKey'], WebEnv=post_result['WebEnv'])
    summary_result=Entrez.read(summary_handle)
    replace_ids={}
    x=0
    for r in summary_result:
        if r['Status']=='replaced':
            x=x+1
            cur_acc=r['ReplacedBy']
            if '.' in cur_acc: cur_acc=cur_acc.replace(cur_acc[cur_acc.find('.'):],'')
            cur_gi=r['Gi']
            if cur_gi in order_ids:
                results[order_ids[cur_gi]]=cur_acc
        elif r['Status']=='suppressed' and remove_old:
            cur_gi=r['Gi']
            if cur_gi in order_ids:
                results[order_ids[cur_gi]]='-'            
    print "total replaced "+str(x)
    return acc_to_gi(results, rettype)
            
    


def get_ids(parts):
    query_ids=[]
    for idx, p in enumerate(parts):
        tmp_p=p
        gi=None
        acc=None
        match_acc=re.search(r'[a-zA-Z]{2}_[0-9]+',p)
        if match_acc: 
            acc=match_acc.group(0)
            tmp_p=tmp_p.replace(acc,'')
        match_gi=re.search(r'([0-9]{4,})',p)
        if match_gi: gi=match_gi.group(0)
        if acc:
            query_ids.append(acc)
        elif gi:
            query_ids.append(gi)

    return query_ids
                


def update_ids(id_list, rettype='acc'):
    query_ids=acc_to_gi(id_list, 'gi')
    results=get_updated(query_ids, rettype)
    return results


def process_text(email, id_text):
    Entrez.email=email
    parts=re.findall(r"[\w']+", id_text)
    id_list=get_ids(parts)
    gene_ids=get_gene_ids(id_list)
    result =""
    for i in gene_ids: result=result+"\t".join(i)+"\n"
    return result

