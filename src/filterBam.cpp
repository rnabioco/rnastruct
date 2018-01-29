#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include <algorithm>

using namespace std ;

void usage() {
  std::cerr << "samtools view -h alignments.bam "
            << "| filterBam -d 5 " << endl
            << "to trim alignments at deletions >= 5 in length "
            << endl;
}

class Cigar {
  public:
    string str ;
    vector<char> cig_ops ;
    vector<int> cig_lens ;
    
    Cigar(){ str = string("") ;}
    Cigar(string x) { str = x ;}
    Cigar(vector<int> x, vector<char> y) {
      cig_lens = x ;
      cig_ops = y ;
      str = this->cigvec_to_str() ;
    }
    
    bool is_aligned() {
      if (str == "*"){
        return false ;
      } else {
        return true ;
      }
    }
    
    void parse_cigar(){
      size_t n;
      char op;
      
      istringstream iss(str) ;
      while(iss >> n >> op){
        cig_lens.push_back(n) ;
        cig_ops.push_back(op) ;
      }
      
      assert(cig_ops.size() == cig_lens.size()) ; 
    }
    
    string cigvec_to_str(){
      size_t n = cig_ops.size() ;
      str = "" ;
      for(int i=0; i<n; ++i){
        str += to_string(cig_lens[i]) ;
        str += cig_ops[i] ;
      }
      return str ; 
    }
    
    int cigar_left_pos(int cigar_op_idx){
      //return distance consumed left of cigar operation
      int consumed_len(0) ;
      for (int i = 0; i < cigar_op_idx; ++i){
        consumed_len += cig_lens[i] ;
      }
      return consumed_len ;
    }
};

class Samline {
  public:
    string readid, contig, cig, rnext, seq, qual, aux;   
    int flag, pos, mapq, pnext, tlen;
}; 

// writer for sam output
std::ostream& operator<< (std::ostream &out, Samline const& data) {
  out << data.readid << '\t' ;
  out << data.flag << '\t' ;
  out << data.contig << '\t' ;
  out << data.pos << '\t' ;
  out << data.mapq << '\t' ;
  out << data.cig << '\t' ;
  out << data.rnext << '\t' ;
  out << data.pnext << '\t' ;
  out << data.tlen << '\t' ;
  out << data.seq << '\t' ;
  out << data.qual << '\t' ;
  out << data.aux << '\t' ; 
  return out ;
}

void modify_sequence(Cigar &cig, 
                     string &oldseq, 
                     string &oldqual,
                     string &newseq, 
                     string &newqual)
{
  size_t n_ops = cig.cig_ops.size() ;
  size_t idx = 0;
  
  for (int i = 0; i < n_ops; ++i){
    char op = cig.cig_ops[i] ;
    int n = cig.cig_lens[i] ;
    switch (op) {
    case 'M':
      newseq += oldseq.substr(idx, n);
      newqual += oldqual.substr(idx, n);
      idx += n;
      break;
    case 'I':
      newseq += oldseq.substr(idx, n);
      newqual += oldqual.substr(idx, n);
      idx += n;
      break;
    case 'D':
      break;
    case 'S':
      newseq += oldseq.substr(idx, n);
      newqual += oldqual.substr(idx, n);
      idx += n;
      break;
    case 'H':
      break;
    case 'P':
      break;
    case '=':
      newseq += oldseq.substr(i, n);
      newqual += oldqual.substr(i, n);
      idx += n;
      break;
    case 'X':
      newseq += oldseq.substr(i, n);
      newqual += oldqual.substr(i, n);
      idx += n;
      break;
    }
  }
}

int find_dels(Cigar cig, int max_del = 4){
  // return first position of cigar deletion operation in cigar vector that exceeds max del
  int idx(0) ;
  size_t n_ops = cig.cig_ops.size() ;
  while(idx < n_ops){
    if(cig.cig_ops[idx] == 'D' && cig.cig_lens[idx] >= max_del) {
      return idx ;
    }
    ++idx ;
  }
  // if no dels found return -1
  return -1 ;
}

bool check_cigar(string cigar_str, int max_del, Cigar &new_cigar){
  bool mod_seq = false; 
  Cigar cig(cigar_str) ;
  
  if (!cig.is_aligned()) {
    return mod_seq;
  }
  cig.parse_cigar() ;
  int del_op = find_dels(cig, max_del) ;

  if (del_op >= 0){
    mod_seq = true ;
    vector<int> c_lens(cig.cig_lens.begin(), cig.cig_lens.begin() + del_op) ;
    vector<char> c_ops(cig.cig_ops.begin(), cig.cig_ops.begin() + del_op) ;
    new_cigar.cig_lens = c_lens ;
    new_cigar.cig_ops = c_ops ;
    new_cigar.cigvec_to_str() ;
  }
  
  return mod_seq ;
}

void split_sam(string line, int max_del, bool clip_read = true){
  // by default clip sequences at deletions for now
  
  Samline sam ;
  stringstream ss(line);
  string flag, pos, mapq, pnext, tlen ;
  getline(ss, sam.readid, '\t');

  getline(ss, flag, '\t');
  sam.flag = stoi(flag) ;
  
  getline(ss, sam.contig, '\t');
  
  getline(ss, pos, '\t');
  sam.pos = stoi(pos) ;
  
  getline(ss, mapq, '\t');
  sam.mapq = stoi(mapq) ;
  
  getline(ss, sam.cig, '\t') ;

  getline(ss, sam.rnext, '\t') ;
  getline(ss, pnext, '\t') ;
  sam.pnext = stoi(pnext) ;
    
  getline(ss, tlen, '\t') ;
  sam.tlen = stoi(tlen) ;
  
  getline(ss, sam.seq, '\t') ;
  getline(ss, sam.qual, '\t') ;
  getline(ss, sam.aux) ;
  
  bool mod_seq ; 
  Cigar new_cigar ;
  mod_seq = check_cigar(sam.cig, max_del, new_cigar) ; 
  
  if(mod_seq){
    if(clip_read){
      string outseq ;
      string outqual ;
      modify_sequence(new_cigar, sam.seq, sam.qual, outseq, outqual) ; 
      sam.seq = outseq ;
      sam.qual = outqual ;
      sam.cig = new_cigar.str ; 
      // report record
      cout << sam << endl ; 
    } else {
      // drop record
      ;
    }
  } else {
    // report unmodified record
    cout << sam << endl ; 
  }
}
  
int main(int argc, char* argv[]){
  
  int del_length = 5 ;
  
  if (argc > 1) {
    std::string arg = argv[1];
    if ((arg == "-h") || (arg == "--help")) {
      usage();
      return 1;
    } else if (arg == "-d") {
      if (argc == 3) { 
        del_length = stoi(argv[2]); 
      } else {
        cerr << "-d requires one argument" << endl;
        usage() ;
        return 1;
      }
    } else {
      usage() ;
      return 1 ;
    }
  } else {
    cerr << "-d option is required" << endl ; 
    usage() ;
    return 1 ;
  }
  
  string line ;
  getline(cin, line) ;
  while(cin) {
    try {
      // skip header
      string str_start = line.substr(0, 1) ;
      if (str_start == "@"){
        cout << line << endl ;
      } else {
        split_sam(line, del_length) ;
      }
      
    } catch(const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
      break;
    }
    getline(cin, line);
  }
}
