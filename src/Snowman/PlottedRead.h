#ifndef SNOWMAN_PLOTTED_READ_H__
#define SNOWMAN_PLOTTED_READ_H__


// hold a single read to be plotted in ASCII plot
struct PlottedRead {

  int pos;
  std::string seq;
  std::string info;

  bool operator<(const PlottedRead& pr) const {
    return (pos < pr.pos);
  }

};

typedef std::vector<PlottedRead> PlottedReadVector;

struct PlottedReadLine {

  std::vector<PlottedRead*> read_vec;
  int available = 0;
  int contig_len = 0;

  void addRead(PlottedRead *r) {
    read_vec.push_back(r);
    available = r->pos + r->seq.length() + 5;
  }

  bool readFits(PlottedRead &r) {
    return (r.pos >= available);
  }

  friend std::ostream& operator<<(std::ostream& out, const PlottedReadLine &r) {
    int last_loc = 0;
    for (auto& i : r.read_vec) {
      assert(i->pos - last_loc >= 0);
      out << std::string(i->pos - last_loc, ' ') << i->seq;
      last_loc = i->pos + i->seq.length();
    }
    int name_buff = r.contig_len - last_loc;
    assert(name_buff < 1e6);
    out << std::string(std::max(name_buff, 5), ' ');
    for (auto& i : r.read_vec) { // add the data
      out << i->info << ",";
    }
    return out;
  }

};

typedef std::vector<PlottedReadLine> PlottedReadLineVector;

#endif
