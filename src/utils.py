import pysam
import shutil
import subprocess
import gzip
import binascii
import os
import sys

def cleanup(tmp_dir, delete = True):
    
    if delete:
      print("removing temp directory: {}".format(tmp_dir))
      shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)     

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return shutil.which(name) is not None



def to_fasta(df, outname, nuc_col_idx = 3, header_str = "1",
        other_data_idx = None):
    """
    write data as fasta format
    Returns: None
    """

    with open(outname, "w") as fout:
        hdr = ">" + header_str
        seq = "".join(df.iloc[:, nuc_col_idx].values)
        
        fout.write(hdr + "\n")
        fout.write(seq + "\n")
        
        if other_data_idx is not None:
            aux_dat = df.iloc[:, other_data_idx].values
            aux_dat = [str(x) for x in aux_dat]
            aux_dat = "".join(aux_dat)
            fout.write(aux_dat + "\n")

def retrieve_contigs(bam, require_reads = True):
    """ return list of contigs, optionally those with reads based on
    index"""
    f = pysam.AlignmentFile(bam, 'rb')

    if require_reads: 
        contigs = []
        for contig in f.get_index_statistics():
            if contig[1] > 0:
                contigs.append(contig[0])  
    else:
        contigs = list(f.references)
    
    f.close()

    return contigs


def gz_is_empty(fname):
    ''' Test if gzip file fname is empty
        Return True if the uncompressed data in fname has zero length
        or if fname itself has zero length
        Raises OSError if fname has non-zero length and is not a gzip file
        https://stackoverflow.com/questions/37874936/how-to-check-empty-gzip-file-in-python
    '''
    if not os.path.isfile(fname) :
      return True
    
    if fname.endswith(".gz"):
        with gzip.open(fname, 'rb') as f:
            data = f.read(1)
        return len(data) == 0
    else: 
       with open(fname, 'r') as f:
            data = f.read(1)
       return len(data) == 0

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def setdiff(lst_a, lst_b):
    sety = set(lst_b)
    return [x for x in lst_a if x not in sety]


def ungzip(in_fn, out_fn):

    with gzip.open(in_fn, 'rb') as f_in:
      with open(out_fn, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    
    os.unlink(in_fn)

def bgzip(in_fn):
    """
    convert file to bgzipped format
    """
    if is_gz_file(in_fn):
      tmp_out_fn = in_fn.replace(".gz", "")
      out_fn = in_fn.replace(".gz", ".bgz")
      ungzip(in_fn, tmp_out_fn)
    else:
      tmp_out_fn = in_fn
      out_fn = in_fn + ".bgz"
    
    pysam.tabix_compress(tmp_out_fn, out_fn, force = True)
    os.unlink(tmp_out_fn)
    
    return out_fn

def unix_sort(out_tmp_ptable, tmp_dir, threads, memory = "8G", verbose = False):


    sort_test = subprocess.run(["sort", "--parallel=2"],
              input = "hello world",
              encoding = 'ascii')
    if sort_test.returncode == 0:
      parallel_sort_available = True
    else:
      parallel_sort_available = False
      if verbose:
        print("Parallel sort not available (man sort), using 1 thread",
               sys.stderr)

    outfile = os.path.join(tmp_dir, "tmp_sorted_pileup_table.tsv")
    outfn = open(outfile, 'w')
    sort_cmd = "gunzip -c " + out_tmp_ptable + \
               " | sort -S " + memory
    if parallel_sort_available:
        sort_cmd = sort_cmd + " --parallel=" + str(threads)

    sort_cmd = sort_cmd + " -k1,1 -k2,2n -k3,3 " + \
               " - "

    if verbose:
        print("sort command is:\n" + sort_cmd, file = sys.stderr)

    sort_run = subprocess.run(sort_cmd, shell=True,
            stderr = sys.stderr,
            stdout = outfn)

    outfn.close()
    return outfile
           
def external_sort(out_tmp_ptable, tmp_dir, nlines, verbose = False):

    if verbose:
        print("started sorting anti-sense and sense pileup tables",
              file = sys.stderr)
    
    nlines_per_chunk, n_chunks = compute_chunks(nlines, 5000000)    
    
    if verbose:
        print("chunking pileup tables into {} lines in {} files".format(nlines_per_chunk, n_chunks),
              file = sys.stderr)
    
    ## chunk into sep. files and sort in memory
    chunk_names = []
    with gzip.open(out_tmp_ptable, 'rt') as input_file:
        header = input_file.readline() 
        for chunk_number in itertools.count(1):
            # read in next chunk of lines and sort them
            lines = itertools.islice(input_file, nlines_per_chunk)
            formatted_lines = []
            for x in lines:
                vals = x.split("\t")
                vals[1] = int(vals[1]) # start
                formatted_lines.append(vals)

            sorted_chunk = sorted(formatted_lines, key = itemgetter(0,1,2))
            if not sorted_chunk:
                # end of input
                break
             
            chunk_name = os.path.join(tmp_dir, 'chunk_{}.chk'.format(chunk_number))
            chunk_names.append(chunk_name)
            with gzip.open(chunk_name, 'wt', compresslevel = 3) as chunk_file:
                for line in sorted_chunk:
                  line[1] = str(line[1])
                  chunk_file.write("\t".join(line))
    
    if verbose:
        print("merging sorted tables",
              file = sys.stderr)
              
    with ExitStack() as stack, gzip.open(out_tmp_ptable, 'wt', compresslevel = 6) as output_file:
        files = [stack.enter_context(gzip.open(chunk, 'rt')) for chunk in chunk_names]
        output_file.write(header)
        output_file.writelines(heapq.merge(*files, key = lambda x: (x.split("\t")[0], 
                                                                    x.split("\t")[1],
                                                                    x.split("\t")[2])))
    for i in chunk_names:
      os.unlink(i)
