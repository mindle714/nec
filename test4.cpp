#include <iostream>
#include <fst/matcher.h>
#include <fst/vector-fst.h>
#include <fst/symbol-table.h>
#include <vector>
#include <string>
#include <fst/determinize.h>
#include <fst/intersect.h>
#include <fst/shortest-path.h>

using SM = fst::SigmaMatcher<
  fst::HashMatcher<fst::Fst<fst::StdArc>>>;

// TODO no consecutive *(play ** by *)
// TODO to be precise, +, not *
int main() {
  std::string pat = "play * by *";

  fst::SymbolTable st;

  st.AddSymbol("_epsilon_");
  for (const auto& c : pat)
    st.AddSymbol(std::string(1, c));
  int sigma = st.Find("*");

  fst::StdVectorFst _fst;

  _fst.AddState();
  _fst.SetStart(0);

  int ps = 0;
  for (int i = 0; i < pat.size(); ++i) {
    if (pat[i] == '*') {
      _fst.AddArc(ps, fst::StdArc(sigma, sigma, 0, ps));
      continue;
    }
    int lab = st.Find(std::string(1, pat[i]));
    int s = _fst.AddState();
    _fst.AddArc(ps, fst::StdArc(lab, lab, 1, s));
    ps = s;
  }

  _fst.SetFinal(ps, 0);

  _fst.Write("one.fst");

  fst::StdVectorFst _fst2;

  std::string input = "play over the rainbow by gundam z";
  for (const auto& c : input)
    st.AddSymbol(std::string(1, c));

  _fst2.AddState();
  _fst2.SetStart(0);

  ps = 0;
  for (int i = 0; i < input.size(); ++i) {
    int lab = st.Find(std::string(1, input[i]));
    int s = _fst2.AddState();
    _fst2.AddArc(ps, fst::StdArc(lab, lab, 0, s));
    ps = s;
  }

  _fst2.SetFinal(ps, 0);

  _fst2.Write("two.fst");

  SM *_mat = new SM(_fst, fst::MATCH_OUTPUT, sigma);
  SM *_mat2 = new SM(_fst2, fst::MATCH_INPUT, sigma);

  fst::CacheOptions copts;
  fst::IntersectFstOptions<fst::StdArc, SM> 
    iopts(copts, _mat, _mat2);

  fst::StdVectorFst _fst3;
  _fst3 = fst::IntersectFst<fst::StdArc>(_fst, _fst2, iopts);
  _fst3.Write("three.fst");

  fst::StdVectorFst _fst4;
  fst::ShortestPath<fst::StdArc>(_fst3, &_fst4);
  _fst4.Write("four.fst");

  if (_fst4.NumStates() == 0) {
    std::cout << "no patterns matched" << std::endl;
    return 0;
  }

  int ss = _fst4.Start();
  while (_fst4.NumArcs(ss) > 0) {
    fst::ArcIterator<fst::StdVectorFst> aiter(_fst4, ss);
    if (aiter.Value().weight == 0) {
      std::cout << st.Find(aiter.Value().ilabel) << std::endl;
    }
    ss = aiter.Value().nextstate;
  }
}
