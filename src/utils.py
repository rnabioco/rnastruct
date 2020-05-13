import pysam
import shutil
import subprocess
import gzip
import binascii
import os
import sys
import pdb
from yaml import load

def is_complex_indel(ref, alt):
    """
    compare two indels and determine if indel is a complex mutation
    https://genome.sph.umich.edu/wiki/Variant_classification
    """ 

    is_complex = False
    rl = len(ref)
    al = len(alt)
    dl = al - rl
    if dl == 0:
        if(rl> 1):
            raise TypeError("MVPs are not supported")
        else:
            raise TypeError("{} is not an indel of {}".format(alt, ref))
     
    # trim variants
    while rl != 1 and al != 1 :
        if ref[rl-1] == alt[al-1]:
            rl -= 1
            al -= 1
        else:
            break

    si = 0

    while rl != 1 and al != 1:
        if ref[si] == alt[si]:
            rl -= 1
            al -= 1
            si += 1
        else:
            break
    ref = ref[si:rl + si]
    alt = alt[si:al + si]
    
    rl = len(ref)
    al = len(alt)

    if rl == 1:
      if ref[0] != alt[0] and  ref[0] != alt[-1]:
        is_complex = True
    elif al == 1:
      if alt[0] != ref[0] and  alt[0] != ref[-1]:
        is_complex = True
    else:
      if dl > 0:
        # insertion
        for i in range(rl):
            if  ref[i] != alt[i]:
                is_complex = True
                break
      else:
          # deletion
        for i in range(al):
            if  ref[i] != alt[i]:
                is_complex = True
                break

    # need to return trim indicies for 5' or 3'
    return is_complex
            
    
def conv_args(pileup_args):
    res = [" "]
    for k,v in pileup_args.items():
        if k == "max_depth":
            res.append("-d")
            res.append(v)
            res.append("-L")
            res.append(v)
        elif k == "ignore_overlaps":
            if v:
                res.append("-x")
        elif k == "min_base_quality":
            res.append("-Q")
            res.append(v)
        elif k == "min_mapping_quality":
            res.append("-q")
            res.append(v)
        elif k == "compute_baq":
            if not v:
                res.append("-B")
        elif k == "ignore_orphans":
            if not v:
                res.append("-A")
        elif k == "redo_baq":
            if v:
                res.append("-E")
        elif k == "adjust_capq_threshold":
            res.append("-C")
            res.append(v)
        else:
            sys.exit(k + " not implemented for now")
    res.append(" ")
    res = [str(x) for x in res]
    return " ".join(res)


def get_pileup_args(fn = None):
  
  if fn is not None:
    config_fn = fn
  else:
    config_fn = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "pileup.yaml")
  with open(config_fn) as f:
      config_data = f.read()
  pileup_args = load(config_data)
  
  return pileup_args
  
pileup_args = get_pileup_args()

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


def ungzip(in_fn, out_fn, remove = True):
    
    with gzip.open(in_fn, 'rb') as f_in:
      with open(out_fn, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    if remove:
      os.unlink(in_fn)

def bgzip(in_fn, remove = True):
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
    if remove:
      os.unlink(tmp_out_fn)
    
    return out_fn

def unix_sort(in_fn, out_fn, threads, memory = "8G", verbose = False):


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
    
    sort_cmd = ["sort", "-S", memory]
    if parallel_sort_available:
        sort_cmd += ["--parallel=" + str(threads)]

    sort_cmd += ["-k1,1", "-k2,2n", "-k3,3"]

    out_fh = gzip.open(out_fn, 'wt', compresslevel = 6)
    
    gunzip = subprocess.Popen(['gunzip', '-c', in_fn],
            stdout=subprocess.PIPE)
    sort_output = subprocess.Popen(sort_cmd,
            stdin = gunzip.stdout,
            stdout = subprocess.PIPE)
    gzip_output = subprocess.Popen(["gzip"],
            stdin = sort_output.stdout,
            stdout = out_fh)
    gunzip.stdout.close()
    gzip_output.wait()

    if verbose:
        print("sort command is:\n" + " ".join(sort_output.args), file = sys.stderr)

    out_fh.close()
    return out_fn
           
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
      



def split_bam(bam, outbam, flags, threads = 1, memory = "2G", force = True):
    """
    bam: bamfile
    outbam: name of output bam
    flags: list of flags to filter bam i.e 
         ["-f 128 -F 16", "-f 80"]
         will generate two tmp bams then merge them
    threads: threads to pass to samtools commands
    memory: memory arg for samtools sort
    force: if outbam exists overwrite
    """
     
    idx_threads = min(threads, 8)
    view_threads = threads
    merge_threads = threads

    if os.path.isfile(outbam):
        if not force:
            sys.exit("{} already exists, set force = True to overwrite".format(outbam))

    tmp_bams = []
    
    for idx,flag_cmd in enumerate(flags):
      tmp_bam = outbam + "_" + str(idx) + ".bam"
      
      view_args = [
              "samtools",
              "view",
              "-b",
              "-@",
              str(view_threads),
              *flag_cmd.split(),
              "-o",
              tmp_bam,
              bam
      ]
      call = subprocess.run(view_args,
              stderr = sys.stderr,
              stdout = sys.stdout)
      tmp_bams.append(tmp_bam)

    merge_args = [
            "samtools",
            "merge",
            "-@",
            str(merge_threads), 
            "-f", 
            outbam, 
            *tmp_bams]

    call = subprocess.run(merge_args,
              stderr = sys.stderr,
              stdout = sys.stdout)

    index_args = [
           "samtools",
           "index",
           "-@",
           str(idx_threads), 
           outbam
    ]

    call = subprocess.run(index_args,
              stderr = sys.stderr,
              stdout = sys.stdout)
    
    if call.returncode != 0:
      
      print("error detected in indexing, trying to sort then reindexing", 
              file = sys.stderr)
      tmpbam = outbam + ".tmp.bam"
      sort_args = [
              "samtools",
              "sort", 
              "-@",
              str(threads),
              "-m",
              memory,
              "-o",
              tmpbam,
              outbam
              ]
      call = subprocess.run(sort_args,
              stderr = sys.stderr,
              stdout = sys.stdout)
      
      if call.returncode != 0:
          sys.exit("problem splitting bam, exiting")
      else:
          os.unlink(outbam)
          os.rename(tmpbam, outbam)

          call = subprocess.run(index_args,
                  stderr = sys.stderr,
                  stdout = sys.stdout)

    for tmp in tmp_bams:
      os.unlink(tmp)

    


