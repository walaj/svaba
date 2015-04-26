#include "BWAWrapper.h"
#include <sstream>
#include <regex>

#define MEM_F_SOFTCLIP  0x200

extern "C" {
  #include <string.h>
  #include "bwamem.h"
}

//mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf);

  /*typedef struct {
  bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
  smem_aux_t *a;
  a = calloc(1, sizeof(smem_aux_t));
  a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
  a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
  return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
  free(a->tmpv[0]->a); free(a->tmpv[0]);
  free(a->tmpv[1]->a); free(a->tmpv[1]);
  free(a->mem.a); free(a->mem1.a);
  free(a);
  }*/

typedef struct {
  bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

#include <iostream>

static inline void makebseq1(const std::string &name, const std::string &seq, bseq1_t *s) {

  s->name = strdup(name.c_str());
  s->comment = 0; //ks->comment.l? strdup(ks->comment.s) : 0;
  s->seq = strdup(seq.c_str());
  s->qual = 0; //ks->qual.l? strdup(ks->qual.s) : 0;
  s->l_seq = seq.length();

}

/*void BWAWrapper::alignSequences(const BWAReadVec &seqs, std::unique_ptr<bwaidx_t>* idx) {

  //new
  for (int i = 0; i < m_blen; i++) {
    mem_alnreg_v reg = mem_align1(memopt, (*idx)->bwt, (*idx)->bns, (*idx)->pac, seqs[i].second.length(), seqs[i].second.c_str());
    mem_aln_t aln = mem_reg2aln(memopt, (*idx)->bns, (*idx)->pac, seqs[i].second.length(), seqs[i].second.c_str(), reg.a);
  }

  }*/

void BWAWrapper::addSequences(const BWAReadVec &seqs, bwaidx_t* idx, SamRecordVec &sam) {

  int m = 0; 

  for (BWAReadVec::const_iterator it = seqs.begin(); it != seqs.end(); it++) {

    // see bseq_read in bwa
    if (m_blen >= m) {
      m = m? m<<1 : 256;
      m_bseqs = (bseq1_t*)realloc(m_bseqs, m * sizeof(bseq1_t));
    }

    makebseq1(it->first, it->second, &m_bseqs[m_blen]);
    m_bseqs[m_blen].id = m_blen;
    m_blen++;

  }

  mem_opt_t * memopt = mem_opt_init();
  memopt->flag |= MEM_F_SOFTCLIP;
  int n_processed;
  mem_pestat_t *fake_mem_pestat_t;

  /*
  //new
  for (int i = 0; i < m_blen; i++) {
    mem_alnreg_v reg = mem_align1(memopt, (*idx)->bwt, (*idx)->bns, (*idx)->pac, seqs[i].second.length(), seqs[i].second.c_str());
    mem_aln_t aln = mem_reg2aln(memopt, (*idx)->bns, (*idx)->pac, seqs[i].second.length(), seqs[i].second.c_str(), reg.a);
    
    //void* aux =  malloc(sizeof(smem_aux_t));
    //smem_aux_t* aux = malloc(sizeof(smem_aux_t));
    //aux = smem_aux_init();
    //char cseq[seqs[i].second.size()+1];//as 1 char space for null is also required
    //strcpy(cseq, seqs[i].second.c_str());
    //mem_align1_core(memopt, (*idx)->bwt, (*idx)->bns, (*idx)->pac, i.second.length(), cseq, aux);
    //mem_alnreg_v reg = mem_align1_core(memopt, (*idx)->bwt, (*idx)->bns, (*idx)->pac, seqs[i].second.length(), cseq, 0);
    //mem_mark_primary_se(memopt, reg.n, reg.a, i);
    //mem_reg2sam(memopt, (*idx)->bns, (*idx)->pac, cseq, &reg, 0, 0);
    //free(reg.a);
    //smem_aux_destroy(aux);
    //free(aux);
  }*/

  mem_process_seqs(memopt, (idx)->bwt, (idx)->bns, (idx)->pac, n_processed, m_blen, m_bseqs, fake_mem_pestat_t);

  for (int i = 0; i < m_blen; i++) {
    //bam1_b * b = malloc(sizeof(bam1_t));
    //sam_parse1(kstring_t *s, bam_hdr_t *h, b)

    // test convert to bam1_t
    //bam1_t *b = bam_init1();
    //b->core

    sam.push_back(SamRecord(std::string(m_bseqs[i].sam)));
  }

  // free all of the char* fields of bseq1_t
  for (int i = 0; i < m_blen; i++) {
    free(m_bseqs[i].seq);
    free(m_bseqs[i].name);
    free(m_bseqs[i].sam);
    free(m_bseqs[i].qual);
    free(m_bseqs[i].comment);;
    m_bseqs[i].seq = NULL;
    m_bseqs[i].name = NULL;
    m_bseqs[i].sam = NULL;
    m_bseqs[i].qual = NULL;
    m_bseqs[i].comment = NULL;

   }

  // free the array and the dummy vars passed to mem_process_seqs
  free (m_bseqs);
  m_bseqs = NULL;
  free (fake_mem_pestat_t);
  free (memopt);

}

/*void BWAWrapper::memProcess(bwaidx_t *idx) {

  int i;
  std::cout << "m_blen: " << m_blen << std::endl;
  for (i = 0; i < m_blen; i++)
    std::cout << std::string(m_bseqs[i].seq) << std::endl;

  mem_opt_t * memopt = mem_opt_init();
  int n_processed;
  mem_pestat_t *fake_mem_pestat_t;

  mem_process_seqs(memopt, idx->bwt, idx->bns, idx->pac, n_processed, m_blen, m_bseqs, fake_mem_pestat_t);
  
  for (i = 0; i < m_blen; i++)
    SamRecord(std::string(m_bseqs[i].sam));

    }*/

SamRecord::SamRecord(std::string sam) {

  record = sam;

  /*
  std::istringstream iss(sam);
  std::string val, line;

  int i = 0;

  while (getline(iss, line, '\n')) {
    std::istringstream issv(line);
    while(getline(issv, val, '\t')) {
      switch(i) {
      case 0 : name = val; break;
      case 1 : try { flag = std::stoi(val); } catch (...) { std::cerr << "Flag Fail on val: " << val << std::endl; } break;
      case 2 : try { chr  = (val == "*") ? -1 : std::stoi(val);} catch (...) { chr = -1; } break;
      case 3 : pos =  std::stoi(val); break;
      case 4 : mapq = std::stoi(val); mapq_vec.push_back(mapq); break;
      case 5 : cigar = val; break;
	//case 7 : 
      case 9 : seq = val; break;
      }
      i++;
      
      // process the tags
      std::regex reg_xp("^XA:Z:(.*)");
      std::smatch match;
      if (std::regex_search(val, match, reg_xp)) {
	xp = match[1].str();
	parseXPString(xp, mapq_vec, nm_vec);
      }

      // process the tags
      std::regex reg_nm("^NM:[A-Za-z]:(.*)");
      if (std::regex_search(val, match, reg_nm)) {
	nm = std::stoi(match[1].str());
	nm_vec.push_back(nm);
      }

    }
    break; // only do the first record of multi-part
  }
  */
}

std::ostream& operator<<(std::ostream &out, const SamRecord &sam) {

  std::string sep = "\t";
  out << sam.name << sep << sam.flag << sep << sam.chr << sep 
      << sam.pos  << sep << sam.mapq << sep << sam.cigar << sep << sam.seq 
      << sep << "XP:" << sam.xp;
  return out;
}


void SamRecord::parseXPString(std::string xp, std::vector<int> &mapq_vec, std::vector<int> &nm_vec) {

  std::istringstream iss(xp);
  std::string token;
  while(getline(iss, token, ';')) {
    std::string record;
    std::istringstream iss2(token);
    int count = 0;
    while(getline(iss2, record, ',')) {
      if (count==0) {
	//chr.push_back(stoi(record));
      } else if (count==3) {
	mapq_vec.push_back(stoi(record));
      } else if (count==4) {
	nm_vec.push_back(stoi(record));
      }
      count++;
    }
  }
  return;

}

