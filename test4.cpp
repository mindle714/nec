#include <iostream>
#include <fst/matcher.h>
#include <fst/vector-fst.h>
#include <fst/symbol-table.h>
#include <vector>
#include <string>
#include <fst/string-weight.h>
#include <fst/determinize.h>
#include <fst/intersect.h>
#include <fst/union.h>
#include <fst/shortest-path.h>

/*class MyArc : public fst::ArcTpl<fst::TropicalWeight> {
public:
  std::string cls = "";

  MyArc() noexcept(std::is_nothrow_default_constructible<Weight>::value) {}

  template <class T>
  MyArc(Label ilabel, Label olabel, T &&weight, StateId nextstate)
    : fst::ArcTpl<fst::TropicalWeight>(ilabel, olabel, weight, nextstate) {}

  MyArc(Label ilabel, Label olabel, StateId nextstate)
    : fst::ArcTpl<fst::TropicalWeight>(ilabel, olabel, nextstate) {}
};*/

//using SW = fst::StringWeight<std::string>;
//using SW = fst::TropicalWeight;
//using SA = fst::ArcTpl<SW>;
//using SA = fst::StringArc<>;
//using SA = MyArc;
using SA = fst::StdArc;
//using SA = fst::ArcTpl<fst::TropicalWeightTpl<int>>;
using SF = fst::VectorFst<SA>;
using SM = fst::SigmaMatcher<
  fst::HashMatcher<fst::Fst<SA>>>;

const std::string LEFT_NE = "{",
  RIGHT_NE = "}";

SF load_pats(fst::SymbolTable& st,
             fst::SymbolTable& nst,
             int* sigma,
             const std::vector<std::string>& pats) {
  fst::StringArc<> sa(1,1,1,1);

  st.AddSymbol("_epsilon_");
  for (const auto& pat : pats) {
    for (const auto& c : pat)
      st.AddSymbol(std::string(1, c));
  }
  *sigma = st.AddSymbol(LEFT_NE);

  nst.AddSymbol("_non_ne_");

  SF _pfst;

  for (const auto& pat : pats) {
    SF _fst;

    _fst.AddState();
    _fst.SetStart(0);

    int ps = 0;
    for (int i = 0; i < pat.size();) {
      if (std::string(1, pat[i]) == LEFT_NE) {
        std::string cls = "";
        ++i;
        while (i < pat.size() && std::string(1, pat[i]) != RIGHT_NE) {
          cls += std::string(1, pat[i]);
          ++i;
        }
        if (i >= pat.size()) {
          std::cerr << "pattern " << pat << "is malformed" << std::endl;
        }
        ++i;

        int clsi = nst.AddSymbol(cls);

        int anys = _fst.AddState();
        _fst.AddArc(ps, SA(*sigma, *sigma, clsi, anys));
        _fst.AddArc(anys, SA(*sigma, *sigma, clsi, anys));
        ps = anys;
        continue;
      }
      int lab = st.Find(std::string(1, pat[i]));
      int s = _fst.AddState();
      _fst.AddArc(ps, SA(lab, lab, 0, s));
      ps = s;
      ++i;
    }

    _fst.SetFinal(ps, 0);
    fst::Union(&_pfst, _fst);
  }

  _pfst.Write("one.fst");
  return _pfst;
}

void find_nes(const fst::SymbolTable& st_,
              const fst::SymbolTable& nst,
              const SF& _pfst,
              std::string input) {
  fst::SymbolTable st = st_;
  int sigma = st_.Find(LEFT_NE);

  for (const auto& c : input)
    st.AddSymbol(std::string(1, c));

  SF _fst2;
  _fst2.AddState();
  _fst2.SetStart(0);

  int ps = 0;
  for (int i = 0; i < input.size(); ++i) {
    int lab = st.Find(std::string(1, input[i]));
    int s = _fst2.AddState();
    _fst2.AddArc(ps, SA(lab, lab, 0, s));
    ps = s;
  }

  _fst2.SetFinal(ps, 0);

  _fst2.Write("two.fst");

  SM *_mat = new SM(_pfst, fst::MATCH_OUTPUT, sigma);
  SM *_mat2 = new SM(_fst2, fst::MATCH_INPUT, sigma);

  fst::CacheOptions copts;
  fst::IntersectFstOptions<SA, SM> 
    iopts(copts, _mat, _mat2);

  SF _fst3;
  _fst3 = fst::IntersectFst<SA>(_pfst, _fst2, iopts);
  _fst3.Write("three.fst");

  SF _fst4;
  fst::ShortestPath<SA>(_fst3, &_fst4);
  _fst4.Write("four.fst");

  if (_fst4.NumStates() == 0) {
    std::cout << "no patterns matched" << std::endl;
    return;
  }

  int pw = 0;
  std::string ne_str = "";
  int ss = _fst4.Start();
  while (_fst4.NumArcs(ss) > 0) {
    fst::ArcIterator<SF> aiter(_fst4, ss);
    if (aiter.Value().ilabel == 0) {
      ss = aiter.Value().nextstate;
      continue; // skip epsilon
    }

    int w = aiter.Value().weight.Value();
    if (w != pw) {
      if (pw != 0 && !ne_str.empty()) {
        std::cout << nst.Find(pw) << " : " << ne_str << std::endl;
      }
      ne_str.clear();
    }
    ne_str += st.Find(aiter.Value().ilabel);
    pw = w;
    ss = aiter.Value().nextstate;
  }
  if (pw != 0 && !ne_str.empty()) {
    std::cout << nst.Find(pw) << " : " << ne_str << std::endl;
  }

}

// TODO no consecutive *(play ** by *)
// TODO no nested ne
int main() {
  std::vector<std::string> pats = {
    "play {song} by {singer}",
    "I'd like to listen {song} from {singer}",
    "can I here {singer}'s {song} or {song}?",
    "play {song} from {singer}"
  };

  fst::SymbolTable st, nst;
  int sigma;

  SF _pfst = load_pats(st, nst, &sigma, pats);

  find_nes(st, nst, _pfst,
      "play gundam z by macross delta");
  find_nes(st, nst, _pfst,
      "I'd like to listen from the hills from the hills");
  find_nes(st, nst, _pfst,
      "can I here Votoms's spirit or what or game or thrown?");
  find_nes(st, nst, _pfst,
      "play gundam z by from the hills");
  find_nes(st, nst, _pfst,
      "play gundam z from by the hills");
}
