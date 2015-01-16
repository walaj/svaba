#include "seqan_tools.h"

/*
void printMultiAlignment(StringVec st) {

  std::cerr << "printMultiAlignment StringVec Length: " << st.size() << std::endl;

  //Align<TSequence> align;
  //resize(rows(align), st.size());

  for (StringVec::iterator it = st.begin(); it != st.end(); it++) {
    for (StringVec::iterator jt = it+1; jt != st.end(); jt++) {

      std::string* p1 = st.data() + (it - st.begin());
      std::string* p2 = st.data() + (jt - st.begin());

      Align<TSequence> align;
      int score = alignSeq2Seq(p1, p2, align);

      unsigned pos_b0 = clippedBeginPosition(row(align, 0));
      unsigned pos_b1 = clippedBeginPosition(row(align, 1));
      unsigned pos_e0 = clippedEndPosition(row(align, 0));
      unsigned pos_e1 = clippedEndPosition(row(align, 1));

      double cutoff = std::max((pos_e0 - pos_b0) * 3.0, 80.0); //3.5 * std::min(it->length(), jt->length());
      //std::cerr << "Clipped Begin: " << pos_b0 << " Clipped End: " << pos_e0 << std::endl;
      //std::cerr << "Clipped Begin: " << pos_b1 << " Clipped End: " << pos_e1 << std::endl;
      //std::cerr << "Score: " << score << " Cutoff " << cutoff << std::endl;

      if (score > cutoff && pos_b0 > pos_b1) {
	//int padlen = it->size() - pos_b0 - jt->size() + 5;
	//padlen = std::max(5, padlen);
	std::cerr << *p1 << std::string(5, ' ') << "Seq" << it-st.begin() << std::endl;
	std::cerr << std::string(pos_b0, ' ') << *p2 << std::string(5, ' ') << "Seq" << jt-st.begin() << std::endl;
	std::cerr << "1------------------------------------------------------------------------------" << std::endl;
      } else if (score > cutoff && pos_b1 >= pos_b0) {
	//int padlen = jt->size() - pos_b1 - it->size() + 5;
	//padlen = std::max(5, padlen);
	std::cerr << *p2 << std::string(5, ' ') << "Seq" << jt-st.begin() << std::endl;
	std::cerr << std::string(pos_b1, ' ') << *p1 << std::string(5, ' ') << "Seq" << it-st.begin() << std::endl;
	std::cerr << "2------------------------------------------------------------------------------" << std::endl;
      } 
    }
  }

}*/

/*
void printMultiAlignment(const ContigVector &contigs) {

  StringVec sv;
  for (ContigVector::const_iterator it = contigs.begin(); it != contigs.end(); it++)
    sv.push_back(it->getSeq());

  printMultiAlignment(sv);

}
*/

// 
int alignSeq2Seq(std::string *cseq, std::string *qseq, TAlign &align) {

    int match = 4;
    int mismatch = -2;
    int gapopen = -4;
    int gapextend = -2;
    int score;
    
    TSequence t_cseq = *cseq;
    TSequence t_qseq = *qseq;

    resize(rows(align), 2);
    assignSource(row(align,0),t_cseq);
    assignSource(row(align,1),t_qseq);
    score = localAlignment(align, Score<int,Simple>(match, mismatch, gapopen, gapextend));
    
    // try the reverse complement
    double cutoff = 20*4;
    if (score < cutoff) {
      std::string cseq_r = *cseq;
      rcomplement(cseq_r);
      t_cseq = cseq_r;

      TAlign align2;
      resize(rows(align2), 2);
      assignSource(row(align2,0),t_cseq);
      assignSource(row(align2,1),t_qseq);
      int score2 = localAlignment(align2, Score<int,Simple>(match, mismatch, gapopen, gapextend));
      
      // if it works, return the complemented cseq and the alignment and score
      if (score2 > cutoff && score2 > score) {
	cseq = &cseq_r;
	align = align2;
	return score;
      }
      
    }
    
    return score;

}

void rcomplement(std::string &a) {

  std::reverse(&a[0], &a[a.size()]);
  std::string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';
}
