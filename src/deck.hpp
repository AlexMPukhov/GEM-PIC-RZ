// ============================================================================
//  deck.hpp -- generic reader for structured ASCII input decks
//
//  Syntax
//  ------
//      # comment (also "!"), anywhere on a line
//      grid {
//        lambda0  = 0.8 um
//        dr_min   = 40 nm           # a unit is required for dimensional values
//        z_points = 0, 50, 5000 um  # lists: comma or blank separated, one unit
//        outdir   = "run 1"         # strings may be quoted
//      }
//      laser { a0 = 2; waist = 8 um }   # ";" separates statements on one line
//
//  * A deck is a sequence of blocks  name { key = value ... }.
//  * Block names and keys are case-insensitive; units are case-sensitive
//    (MeV vs meV).
//  * A block type may occur several times (laser, bunch); blocks of the same
//    type are addressed as laser[0], laser[1], ... in overrides and messages.
//  * The deck reader only knows the syntax. Which blocks and keys exist, and
//    their units, is defined by the consumer (config.cpp). Every key that the
//    consumer never asks for is reported as an error, so typos are caught.
// ============================================================================
#pragma once
#include <stdexcept>
#include <string>
#include <vector>

namespace gem {

struct DeckError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct Entry {
  std::string key;                  // lower case
  std::vector<std::string> tokens;  // value tokens
  int line = 0;
  mutable bool used = false;        // set when the consumer reads the key
};

struct Block {
  std::string name;  // lower case
  int line = 0;      // line of the block header
  int index = 0;     // index among blocks with the same name
  int count = 1;     // number of blocks with this name in the deck
  std::vector<Entry> entries;

  const Entry* find(const std::string& key) const;
  Entry* find(const std::string& key);
  std::string label() const;  // "grid" or "laser[1]" (index shown if count > 1)
};

class Deck {
 public:
  std::string source = "<string>";  // file name, for messages
  std::vector<Block> blocks;

  static Deck parseFile(const std::string& path);
  static Deck parseString(const std::string& text, const std::string& source = "<string>");

  // all blocks with the given name, in order of appearance
  std::vector<const Block*> all(const std::string& name) const;
  // the single block with this name, or nullptr; error if it occurs twice
  const Block* single(const std::string& name) const;

  // error for any block whose name is not in `known`
  void checkBlockNames(const std::vector<std::string>& known) const;
  // error listing every key that was never read by the consumer
  void checkAllUsed() const;

  // command-line override "block.key=value" or "block[i].key=value";
  // replaces the value, or adds the key if it is not in the block; a block
  // type that does not occur in the deck at all is created
  void applyOverride(const std::string& spec);

  // the deck in canonical form (after overrides), for the run directory
  std::string dump() const;

  std::string where(int line) const { return source + ":" + std::to_string(line); }

 private:
  void finalize();  // set Block::index / Block::count
};

}  // namespace gem
