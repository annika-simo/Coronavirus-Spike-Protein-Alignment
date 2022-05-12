#include "unit_test.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <vector>

/* References:
 * https://bioboot.github.io/bimm143_W20/class-material/nw/
 * https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
 * https://www.ncbi.nlm.nih.gov/datasets/coronavirus/proteins/
 */

using direction = enum : uint8_t { diagonal = 0, horizontal = 1, vertical = 2 };
using NWAlignmentDifference = std::vector<std::vector<int>>;
using NWAlignmentDirection = std::vector<std::vector<direction>>;
using ProteinSequence = std::vector<char>;
using Comparison = std::vector<std::tuple<int, char, char>>;

/*
 * get_id_and_sequence
 *
 * If it is possible to read another id/protein sequence from file, id is 
 * that sequence id, sequence is that sequence and the function returns true. 
 * Otherwise, the function returns false and id and sequence are unchanged.
 *
 * input: std::ifstream, named file, by reference; std::string, named id, by 
 *        reference; std::string, named sequence, by reference
 * output: bool
 * 
 */

bool get_id_and_sequence(std::ifstream &file, std::string &id, std::string &sequence) {
    if (file.eof() || !file.is_open()) {
        return false;
    } else {
        std::getline(file, id);
        std::getline(file, sequence);
        return true;
    }
}

/*
 * open_file
 *
 * If a std::ifstream can be opened for the file named filename, then file is that 
 * std::ifstream and the function returns true. Otherwise, the function returns 
 * false and file is unchanged.
 *
 * input: std::string, named filename; std::ifstream, named file, by reference
 * output: bool
 *
 */

bool open_file(std::string filename, std::ifstream &file) { 
    file.open(filename);
    if (!file.is_open()) {
        return false;
    }
    else {
        return true; 
    }
}

/*
 * string_to_protein_sequence
 *
 * A ProteinSequence that represents that sequence contained in protein_str.
 *
 * input: string named protein_str
 * output: ProteinSequence
 *
 */

ProteinSequence string_to_protein_sequence(const std::string &protein_str) {
    ProteinSequence endSequence;
    for (int i = 0; i < protein_str.length(); ++i) {
        endSequence.push_back(protein_str.at(i));
    }

    return endSequence;
}

/*
 * above_left
 *
 * If (row-1)and (col-1) are greater-than-or-equal-to zero, return the value in the 
 * cell at index (row-1), (col-1).  Otherwise, if col equals 0, return -1 * row. 
 * Otherwise, return -1 * col. above_left will never be given a row/col pair that is 
 * less than 0 or greater than the bounds of difference_vector.
 *
 * input: NWAlignment Difference named difference_vector; row index, named row; column 
          index, named col
 * output: int
 *
 */

int above_left(const NWAlignmentDifference &difference_vector, int row, int col) {
    if ((row - 1) >= 0 && (col - 1) >= 0) {
        return difference_vector[row - 1][col - 1];
    }
    else if (col == 0) {
        return (-1 * row);
    }
    else {
        return (-1 * col);
    }
}

/*
 * above
 *
 * If (row-1) is greater-than-or-equal-to zero, return the value in the cell at index 
 * (row-1), col. Otherwise, return -1 * (col + 1). above will never be given a row/col 
 * pair that is less than 0 or greater than the bounds of difference_vector.	
 *
 * input: NWAlignment Difference named difference_vector; row index, named row; column
          index, named col
 * output: int
 *
 */

int above(const NWAlignmentDifference &difference_vector, int row, int col) {
    if ((row - 1) >= 0) {
        return difference_vector[row - 1][col];
    }
    else {
        return (-1 * (col + 1));
    }
}

/*
 * left
 *
 * If (col-1) is greater-than-or-equal-to zero, return the value in the cell at index row, 
 * (col-1) . Otherwise, return -1 * (row + 1). left will never be given a row/col pair that 
 * is less than 0 or greater than the bounds of difference_vector.	
 *
 * input: NWAlignment Difference named difference_vector; row index, named row; column
          index, named col
 * output: int
 *
 */

int left(const NWAlignmentDifference &difference_vector, int row, int col) {
    if ((col - 1) >= 0) {
        return difference_vector[row][col - 1];
    }
    else {
        return (-1 * (row + 1));
    }
}

void print_comparison(const Comparison &comparison,
                      const std::string &comparison_id, bool verbose = false) {

  std::cout << comparison_id << ": " << comparison.size();

  if (verbose) {
    for (auto i : comparison) {
      auto index = std::get<0>(i);
      auto first = std::get<1>(i);
      auto second = std::get<2>(i);
      std::cout << "(@" << index << ": " << first << ", " << second << "), ";
    }
  }
  std::cout << "\n";
}

template <typename T>
std::vector<T> reverse_vector(const std::vector<T> &to_reverse) {
  std::vector<T> forward;
  for (auto i : to_reverse) {
    forward.push_back(i);
  }
  return forward;
}

Comparison create_comparison(const ProteinSequence &top,
                             const ProteinSequence &bottom,
                             const NWAlignmentDirection &dirs) {
  int row = bottom.size() - 1;
  int col = top.size() - 1;
  Comparison full_result, diff_result;
  while (row >= 0 && col >= 0) {
    char first{0};
    char second{0};
    if (dirs[row][col] == diagonal) {
      first = bottom[row];
      second = top[col];
      row--;
      col--;
    } else if (dirs[row][col] == horizontal) {
      first = '-';
      second = top[col];
      col--;
    } else {
      first = bottom[row];
      second = '-';
      row--;
    }
    full_result.push_back(std::make_tuple(0, first, second));
  }
  full_result = reverse_vector(full_result);

  auto index = 0;
  for (auto i : full_result) {
    auto first = std::get<1>(i);
    auto second = std::get<2>(i);
    if (first != second) {
      diff_result.push_back(std::make_tuple(index, first, second));
    }
    index++;
  }
  return diff_result;
}

NWAlignmentDirection new_direction_vector(int rows, int cols) {
  std::vector<std::vector<direction>> new_direction_vector(
      rows, std::vector<direction>(cols));
  return new_direction_vector;
}

NWAlignmentDifference new_difference_vector(int rows, int cols) {
  NWAlignmentDifference new_diff_vector(rows, std::vector<int>(cols));
  return new_diff_vector;
}

int s(const char a, const char b) {
  if (a == 'X' || b == 'X') {
    return 0;
  }
  if (a == b) {
    return 1;
  }
  return -1;
}

NWAlignmentDirection needleman_wunsch(const ProteinSequence &top, const ProteinSequence &bottom) {
  auto bottom_length = bottom.size();
  auto top_length = top.size();
  auto directions = new_direction_vector(bottom_length, top_length);
  auto differences = new_difference_vector(bottom_length, top_length);
  for (int r = 0; r < bottom_length; r++) {
    for (int c = 0; c < top_length; c++) {
      int t, l, d, match, insert, del;
      t = above(differences, r, c);
      l = left(differences, r, c);
      d = above_left(differences, r, c);

      match = d + s(top[c], bottom[r]);
      del = t - 1;
      insert = l - 1;

      if (match > del && match > insert) {
        differences[r][c] = match;
        directions[r][c] = diagonal;
      } else if (del > insert) {
        differences[r][c] = del;
        directions[r][c] = vertical;
      } else {
        differences[r][c] = insert;
        directions[r][c] = horizontal;
      }
    }
  }
  return directions;
}

void run_unit_tests() {
  std::ifstream test_input_file;

  check_result(true, open_file("testing.fasta", test_input_file));
  check_result(true, test_input_file.is_open());

  std::string id, sequence;

  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id1"});
  check_result(true, sequence == std::string{"sequence1"});
  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id2"});
  check_result(true, sequence == std::string{"sequence2"});

  test_input_file.close();

  NWAlignmentDifference difference = new_difference_vector(2, 2);
  int counter{0};
  for (auto row : difference) {
    for (auto cell : row) {
      cell = counter++;
    }
  }

  check_result(0, above(difference, 1, 0));
  check_result(-2, above(difference, 0, 1));
  check_result(0, left(difference, 0, 1));
  check_result(-6, left(difference, 5, 0));
  check_result(0, above_left(difference, 1, 1));
  check_result(-4, above_left(difference, 0, 4));
  check_result(-3, above_left(difference, 3, 0));
}
 
int main() {
    std::ifstream referenceFile; 
    std::ifstream comparisonFile;
    std::string referenceString;
    std::string id;
    std::string comparisonString;

    bool verbose = false;

    open_file("reference.fasta", referenceFile);
    open_file("comparison.fasta", comparisonFile);

    get_id_and_sequence(referenceFile, id, referenceString);

    ProteinSequence referenceSequence = string_to_protein_sequence(referenceString);

    while (get_id_and_sequence(comparisonFile, id, comparisonString)) {
        if (id == "") {
            break;
        }
        ProteinSequence comparisonSequence = string_to_protein_sequence(comparisonString);
        NWAlignmentDirection direction = needleman_wunsch(referenceSequence, comparisonSequence);
        Comparison comparison = create_comparison(referenceSequence, comparisonSequence, direction);
        print_comparison(comparison, id, verbose);
    }
    // while loop read file - condition return true for id fuction
    // read with fuction
    // return 
    return 0;
}
