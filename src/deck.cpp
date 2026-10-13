// ============================================================================
//  deck.cpp -- tokenizer and parser for the input deck (see deck.hpp)
// ============================================================================
#include "deck.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>

namespace gem {

namespace {

enum class Tok { Word, Str, LBrace, RBrace, Eq, End };

struct Token {
  Tok type;
  std::string text;
  int line;
};

std::string lower(std::string s) {
  for (auto& ch : s) ch = char(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}

bool isSpecial(char ch) {
  return ch == '{' || ch == '}' || ch == '=' || ch == ';' || ch == '#' || ch == '!' || ch == ',' ||
         ch == '"' || ch == '\'';
}

std::vector<Token> tokenize(const std::string& text, const std::string& src) {
  std::vector<Token> out;
  int line = 1;
  size_t i = 0;
  const size_t n = text.size();
  auto push = [&](Tok t, std::string s) { out.push_back({t, std::move(s), line}); };
  while (i < n) {
    const char ch = text[i];
    if (ch == '\n') { push(Tok::End, "\\n"); ++line; ++i; continue; }
    if (ch == ';') { push(Tok::End, ";"); ++i; continue; }
    if (ch == '#' || ch == '!') { while (i < n && text[i] != '\n') ++i; continue; }
    if (std::isspace(static_cast<unsigned char>(ch)) || ch == ',') { ++i; continue; }
    if (ch == '{') { push(Tok::LBrace, "{"); ++i; continue; }
    if (ch == '}') { push(Tok::RBrace, "}"); ++i; continue; }
    if (ch == '=') { push(Tok::Eq, "="); ++i; continue; }
    if (ch == '"' || ch == '\'') {
      size_t j = i + 1;
      while (j < n && text[j] != ch && text[j] != '\n') ++j;
      if (j >= n || text[j] != ch)
        throw DeckError(src + ":" + std::to_string(line) + ": unterminated string");
      push(Tok::Str, text.substr(i + 1, j - i - 1));
      i = j + 1;
      continue;
    }
    size_t j = i;
    while (j < n && !std::isspace(static_cast<unsigned char>(text[j])) && !isSpecial(text[j])) ++j;
    push(Tok::Word, text.substr(i, j - i));
    i = j;
  }
  push(Tok::End, "<eof>");
  return out;
}

}  // namespace

// ---------------------------------------------------------------- Block
const Entry* Block::find(const std::string& key) const {
  for (const auto& e : entries)
    if (e.key == key) return &e;
  return nullptr;
}
Entry* Block::find(const std::string& key) {
  for (auto& e : entries)
    if (e.key == key) return &e;
  return nullptr;
}
std::string Block::label() const {
  return count > 1 ? name + "[" + std::to_string(index) + "]" : name;
}

// ---------------------------------------------------------------- parsing
Deck Deck::parseFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw DeckError("cannot open input deck '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return parseString(ss.str(), path);
}

Deck Deck::parseString(const std::string& text, const std::string& source) {
  Deck d;
  d.source = source;
  const std::vector<Token> T = tokenize(text, source);
  size_t k = 0;
  auto err = [&](int line, const std::string& msg) -> DeckError {
    return DeckError(source + ":" + std::to_string(line) + ": " + msg);
  };
  auto skipEnd = [&] { while (k < T.size() && T[k].type == Tok::End) ++k; };

  while (true) {
    skipEnd();
    if (k >= T.size()) break;
    if (T[k].type != Tok::Word) throw err(T[k].line, "expected a block name, found '" + T[k].text + "'");
    Block b;
    b.name = lower(T[k].text);
    b.line = T[k].line;
    ++k;
    skipEnd();
    if (k >= T.size() || T[k].type != Tok::LBrace)
      throw err(b.line, "expected '{' after block name '" + b.name + "'");
    ++k;
    bool closed = false;
    while (k < T.size()) {
      skipEnd();
      if (k >= T.size()) break;
      if (T[k].type == Tok::RBrace) { ++k; closed = true; break; }
      if (T[k].type != Tok::Word)
        throw err(T[k].line, "expected 'key = value' in block '" + b.name + "', found '" + T[k].text + "'");
      Entry e;
      e.key = lower(T[k].text);
      e.line = T[k].line;
      ++k;
      if (k < T.size() && T[k].type == Tok::LBrace)
        throw err(e.line, "unexpected '{': block '" + b.name + "' (line " + std::to_string(b.line) +
                              ") is not closed with '}' before block '" + e.key + "'");
      if (k >= T.size() || T[k].type != Tok::Eq)
        throw err(e.line, "expected '=' after key '" + e.key + "' in block '" + b.name + "'");
      ++k;
      while (k < T.size() && (T[k].type == Tok::Word || T[k].type == Tok::Str)) {
        e.tokens.push_back(T[k].text);
        ++k;
      }
      if (e.tokens.empty()) throw err(e.line, "missing value for key '" + e.key + "'");
      if (k < T.size() && T[k].type == Tok::Eq)
        throw err(T[k].line, "unexpected '=' in the value of '" + e.key + "' (missing line break or ';'?)");
      if (k < T.size() && T[k].type == Tok::LBrace)
        throw err(T[k].line, "unexpected '{' in block '" + b.name + "' (missing '}' before a new block?)");
      if (const Entry* old = b.find(e.key))
        throw err(e.line, "key '" + e.key + "' given twice in block '" + b.name + "' (first in line " +
                              std::to_string(old->line) + ")");
      b.entries.push_back(std::move(e));
    }
    if (!closed) throw err(b.line, "block '" + b.name + "' is not closed with '}'");
    d.blocks.push_back(std::move(b));
  }
  d.finalize();
  return d;
}

void Deck::finalize() {
  std::map<std::string, int> cnt;
  for (auto& b : blocks) b.index = cnt[b.name]++;
  for (auto& b : blocks) b.count = cnt[b.name];
}

// ---------------------------------------------------------------- queries
std::vector<const Block*> Deck::all(const std::string& name) const {
  std::vector<const Block*> v;
  for (const auto& b : blocks)
    if (b.name == name) v.push_back(&b);
  return v;
}

const Block* Deck::single(const std::string& name) const {
  auto v = all(name);
  if (v.empty()) return nullptr;
  if (v.size() > 1)
    throw DeckError(where(v[1]->line) + ": block '" + name + "' may appear only once (first in line " +
                    std::to_string(v[0]->line) + ")");
  return v[0];
}

void Deck::checkBlockNames(const std::vector<std::string>& known) const {
  for (const auto& b : blocks) {
    if (std::find(known.begin(), known.end(), b.name) == known.end()) {
      std::string list;
      for (const auto& k : known) list += (list.empty() ? "" : ", ") + k;
      throw DeckError(where(b.line) + ": unknown block '" + b.name + "' (known blocks: " + list + ")");
    }
  }
}

void Deck::checkAllUsed() const {
  std::string msg;
  for (const auto& b : blocks)
    for (const auto& e : b.entries)
      if (!e.used) msg += "\n  " + where(e.line) + ": unknown key '" + e.key + "' in block '" + b.label() + "'";
  if (!msg.empty()) throw DeckError("input deck contains unknown keys:" + msg);
}

// ---------------------------------------------------------------- overrides
void Deck::applyOverride(const std::string& spec) {
  const auto eq = spec.find('=');
  const auto dot = spec.find('.');
  if (eq == std::string::npos || dot == std::string::npos || dot > eq)
    throw DeckError("override '" + spec + "' must have the form block.key=value or block[i].key=value");
  std::string bname = lower(spec.substr(0, dot));
  const std::string key = lower(spec.substr(dot + 1, eq - dot - 1));
  const std::string value = spec.substr(eq + 1);
  int idx = 0;
  const auto br = bname.find('[');
  if (br != std::string::npos) {
    const auto close = bname.find(']', br);
    if (close == std::string::npos) throw DeckError("override '" + spec + "': missing ']'");
    idx = std::stoi(bname.substr(br + 1, close - br - 1));
    bname = bname.substr(0, br);
  }
  Block* target = nullptr;
  for (auto& b : blocks)
    if (b.name == bname && b.index == idx) target = &b;
  if (!target && idx == 0 && all(bname).empty()) {
    // the deck has no block of this type: create it (e.g. numerics.seed=3)
    Block nb;
    nb.name = bname;
    nb.line = 0;
    blocks.push_back(nb);
    finalize();
    target = &blocks.back();
  }
  if (!target)
    throw DeckError("override '" + spec + "': the deck has no block " + bname + "[" + std::to_string(idx) + "]");
  // tokenize the value with the deck tokenizer
  Deck tmp = parseString("x { v = " + value + " }", "override '" + spec + "'");
  const Entry& ve = tmp.blocks.at(0).entries.at(0);
  if (Entry* e = target->find(key)) {
    e->tokens = ve.tokens;
  } else {
    Entry e2;
    e2.key = key;
    e2.tokens = ve.tokens;
    e2.line = 0;
    target->entries.push_back(e2);
  }
}

std::string Deck::dump() const {
  std::ostringstream os;
  for (const auto& b : blocks) {
    os << b.name << " {\n";
    size_t w = 0;
    for (const auto& e : b.entries) w = std::max(w, e.key.size());
    for (const auto& e : b.entries) {
      os << "  " << e.key << std::string(w - e.key.size(), ' ') << " =";
      for (const auto& t : e.tokens) {
        const bool quote = t.empty() || std::any_of(t.begin(), t.end(), [](char ch) {
          return std::isspace(static_cast<unsigned char>(ch)) || isSpecial(ch);
        });
        os << ' ' << (quote ? "\"" + t + "\"" : t);
      }
      os << '\n';
    }
    os << "}\n\n";
  }
  return os.str();
}

}  // namespace gem
