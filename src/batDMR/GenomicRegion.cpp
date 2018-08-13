/*
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "GenomicRegion.hpp"

#include <cassert>
#include <fstream>
#include <tr1/unordered_map>
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

using std::string;
using std::vector;
using std::ostringstream;
using std::tr1::unordered_map;
unordered_map<string, chrom_id_type> SimpleGenomicRegion::fw_table_in;
unordered_map<chrom_id_type, string> SimpleGenomicRegion::fw_table_out;

string strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

chrom_id_type SimpleGenomicRegion::assign_chrom(const std::string &c) {
  unordered_map<string, chrom_id_type>::const_iterator
    chr_id(fw_table_in.find(c));
  if (chr_id == fw_table_in.end()) {
    const chrom_id_type r = fw_table_in.size();
    fw_table_in[c] = r;
    fw_table_out[r] = c;
    // return r;
  }
  // else return chr_id->second;
  return fw_table_in[c];
}

bool is_valid_filename(const string name, const string& filename_suffix) {
  const string suffix(name.substr(name.find_last_of(".") + 1));
  return (suffix == filename_suffix);
}
string path_join(const string& a, const string& b) {
  if (b.empty() || b[0] == '/')
    throw DNAmethException("cannot prepend dir to file \"" + b + "\"");
  if (!a.empty() && a[a.length() - 1] == '/')
    return a + b;
  else
    return a + "/" + b;
}
void read_dir(const string& dirname, vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw "could not open directory: " + dirname;

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw "error reading directory: " + dirname;
  if (filenames.empty())
    throw "no valid files found in: " + dirname;
  closedir(dir);
}
void read_dir(const string& dirname, string filename_suffix,
    vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw DNAmethException("could not open directory: " + dirname);

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw DNAmethException("error reading directory: " + dirname);
  if (filenames.empty())
    throw DNAmethException("no valid files found in: " + dirname);
  closedir(dir);
}

char complement(int i) {
  static const int b2c_size = 20;
  static const char b2c[] = {
   //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
    'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
  };
  static const char b2cl[] = {
   //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
    't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
  };
  if (i - 'A' >= 0 && i - 'A' < b2c_size)
    return b2c[i - 'A'];
  else if (i - 'a' >= 0 && i - 'a' < b2c_size)
    return b2cl[i - 'a'];
  else return 'N';
}
std::string strip(const std::string& s) {
  const size_t len = s.length();
  size_t i = 0;
  while (i < len && isspace(s[i])) {
    i++;
  }
  size_t j = len;
  do {
    j--;
  } while (j >= i && isspace(s[j]));
  j++;
  if (i == 0 && j == len)
    return s;
  else return s.substr(i, j - i);
}
inline void revcomp_inplace(std::string& s) {
  std::transform(s.begin(), s.end(), s.begin(), complement);
  std::reverse(s.begin(), s.end());
}

inline void revcomp_inplace(std::string::iterator first, std::string::iterator last) {
  std::transform(first, last, first, complement);
  std::reverse(first, last);
}

string
SimpleGenomicRegion::retrieve_chrom(chrom_id_type i) {
  unordered_map<chrom_id_type, string>::const_iterator chr_name(fw_table_out.find(i));
  assert(chr_name != fw_table_out.end());
  return chr_name->second;
}

std::vector<std::string> split_whitespace_quoted(std::string to_split) {
  static const char *non_word_chars = " \t\"'";
  to_split = strip(to_split);
  
  std::vector<std::string> words;
  size_t start_pos = 0, end_pos = 0;
  const size_t length_of_to_split = to_split.length();
  while (start_pos < length_of_to_split) {
    /** find next position that is not a word character */
    end_pos = to_split.find_first_of(non_word_chars, end_pos);
    if (end_pos > to_split.length()) { /** If we hit the end: done */
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      break;
    }
    /** unescaped, unquoted white space: definitely a word delimiter */
    if (to_split[end_pos] == ' ' || to_split[end_pos] == '\t') { 
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      end_pos = to_split.find_first_not_of(" \t", end_pos);
      start_pos = end_pos;
    }
    /** preserve whatever is being escaped; will become part of the
	current word */
    else if (to_split[end_pos] == '\\')
      end_pos = to_split.find_first_not_of(non_word_chars, end_pos + 2);
    else {
      const std::string current_delim = "\\" + to_split.substr(end_pos, 1);
      do { // slurp doubly- or singly-quoted string
	end_pos = to_split.find_first_of(current_delim, end_pos + 1);
	if (end_pos == std::string::npos) {
	  end_pos = length_of_to_split;
	  break;
	}
	if (to_split[end_pos] == '\\')
	  ++end_pos;
	else break;
      } while (true);
      ++end_pos;
    }
    if (end_pos >= length_of_to_split) {
      words.push_back(to_split.substr(start_pos,
				      end_pos - start_pos));
      break;
    }
  }
  return words;
}

SimpleGenomicRegion::SimpleGenomicRegion(const GenomicRegion &r) :
  chrom(assign_chrom(r.get_chrom())), start(r.get_start()), end(r.get_end()) {}

SimpleGenomicRegion::SimpleGenomicRegion(string string_representation) {
  vector<string> parts = split_whitespace_quoted(string_representation);

  // make sure there is the minimal required info
  if (parts.size() < 3)
    throw GenomicRegionException("Invalid string representation: " +
                                 string_representation);
  // set the chromosome name
  chrom = assign_chrom(parts[0]);

  // set the start position
  const int checkChromStart = atoi(parts[1].c_str());
  if (checkChromStart < 0)
    throw GenomicRegionException("Invalid start: " + parts[1]);
  else start = static_cast<size_t>(checkChromStart);

  // set the end position
  const int checkChromEnd = atoi(parts[2].c_str());
  if (checkChromEnd < 0)
    throw GenomicRegionException("Invalid end: " + parts[2]);
  else end = static_cast<size_t>(checkChromEnd);
}

SimpleGenomicRegion::SimpleGenomicRegion(const char *s, const size_t len) {
  size_t i = 0;

  // the chrom
  while (isspace(s[i]) && i < len) ++i;
  size_t j = i;
  while (!isspace(s[i]) && i < len) ++i;
  chrom = assign_chrom(string(s + j, i - j));

  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  start = atoi(s + j);

  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  end = atoi(s + j);
}

string
SimpleGenomicRegion::tostring() const {
  std::ostringstream s;
  s << retrieve_chrom(chrom) << "\t" << start << "\t" << end;
  return s.str();
}


std::ostream&
operator<<(std::ostream& s, const SimpleGenomicRegion& region) {
  return s << region.tostring();
}

std::istream&
operator>>(std::istream& s, SimpleGenomicRegion& region) {
  string chrom;
  size_t start = 0ul, end = 0ul;
  if (s >> chrom >> start >> end) {
    region = SimpleGenomicRegion(chrom, start, end);
    // else region = SimpleGenomicRegion();
  }
  else s.setstate(std::ios::badbit);

  char c;
  while ((c = s.get()) != '\n' && s);
  if (c != '\n')
    s.setstate(std::ios::badbit);

  if (s.eof())
    s.setstate(std::ios::badbit);

  return s;
}


bool
SimpleGenomicRegion::contains(const SimpleGenomicRegion& other) const {
  return chrom == other.chrom && start <= other.start && other.end <= end;
}


bool
SimpleGenomicRegion::overlaps(const SimpleGenomicRegion& other) const {
  return chrom == other.chrom &&
    ((start < other.end && other.end <= end) ||
     (start <= other.start && other.start < end) ||
     other.contains(*this));
}


size_t
SimpleGenomicRegion::distance(const SimpleGenomicRegion& other) const {
  if (chrom != other.chrom)
    return std::numeric_limits<size_t>::max();
  else if (overlaps(other) || other.overlaps(*this))
    return 0;
  else return (end < other.start) ?
         other.start - end + 1 : start - other.end + 1;
}


bool
SimpleGenomicRegion::operator<(const SimpleGenomicRegion& rhs) const {
  return (get_chrom() < rhs.get_chrom() ||
          (chrom == rhs.chrom &&
           (start < rhs.start ||
            (start == rhs.start && (end < rhs.end)))));
}

bool
SimpleGenomicRegion::less1(const SimpleGenomicRegion& rhs) const {
  return (get_chrom() < rhs.get_chrom() ||
          (chrom == rhs.chrom &&
           (end < rhs.end ||
            (end == rhs.end && start < rhs.start))));
}

bool
SimpleGenomicRegion::operator<=(const SimpleGenomicRegion& rhs) const {
  return !(rhs < *this);
}


bool
SimpleGenomicRegion::operator==(const SimpleGenomicRegion& rhs) const {
  return (chrom == rhs.chrom && start == rhs.start && end == rhs.end);
}


bool
SimpleGenomicRegion::operator!=(const SimpleGenomicRegion& rhs) const {
  return (chrom != rhs.chrom || start != rhs.start || end != rhs.end);
}

#include <iostream>

unordered_map<string, chrom_id_type> GenomicRegion::fw_table_in;
unordered_map<chrom_id_type, string> GenomicRegion::fw_table_out;

chrom_id_type
GenomicRegion::assign_chrom(const std::string &c) {
  unordered_map<string, chrom_id_type>::const_iterator chr_id(fw_table_in.find(c));
  if (chr_id == fw_table_in.end()) {
    const chrom_id_type r = fw_table_in.size();
    fw_table_in[c] = r;
    fw_table_out[r] = c;
    return r;
  }
  else return chr_id->second;
}
using std::cerr;
using std::endl;

string
GenomicRegion::retrieve_chrom(chrom_id_type i) {
  unordered_map<chrom_id_type, string>::const_iterator chr_name(fw_table_out.find(i));
  assert(chr_name != fw_table_out.end());
  return chr_name->second;
}


GenomicRegion::GenomicRegion(string string_representation) : strand('+') {
  vector<string> parts(split_whitespace_quoted(string_representation));

  // make sure there is the minimal required info
  if (parts.size() < 3)
    throw GenomicRegionException("Invalid string representation: " +
                                 string_representation);
  // set the chromosome name
  chrom = assign_chrom(parts[0]);

  // set the start position
  const int checkChromStart = atoi(parts[1].c_str());
  if (checkChromStart < 0)
    throw GenomicRegionException("Invalid start: " + parts[1]);
  else start = static_cast<size_t>(checkChromStart);

  // set the end position
  const int checkChromEnd = atoi(parts[2].c_str());
  if (checkChromEnd < 0)
    throw GenomicRegionException("Invalid end: " + parts[2]);
  else end = static_cast<size_t>(checkChromEnd);

  if (parts.size() > 3)
    name = parts[3];

  if (parts.size() > 4)
    score = atof(parts[4].c_str());

  if (parts.size() > 5)
    strand = parts[5][0];
}


GenomicRegion::GenomicRegion(const char *s, const size_t len) {
  size_t i = 0;

  // the chrom
  while (isspace(s[i]) && i < len) ++i;
  size_t j = i;
  while (!isspace(s[i]) && i < len) ++i;
  chrom = assign_chrom(string(s + j, i - j));

  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  start = atoi(s + j);
  while (!isspace(s[i]) && i < len) ++i;

  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len) ++i;

  // name of the region
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  name = string(s + j, i - j);

  // score of the region (floating point)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len) ++i;
  score = atof(s + j);

  // strand
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len) ++i;
  strand = *(s + j);

  // ADS: This is a hack!!!
  if (strand != '-') strand = '+';
}

string
GenomicRegion::tostring() const {
  std::ostringstream s;
  s << retrieve_chrom(chrom) << "\t" << start << "\t" << end;
  if (!name.empty())
    s << "\t" << name << "\t" << score << "\t" << strand;
  return s.str();
}

std::ostream&
operator<<(std::ostream& s, const GenomicRegion& region) {
  return s << region.tostring();
}

std::istream&
operator>>(std::istream& s, GenomicRegion& region) {
  string chrom, name;
  size_t start = 0ul, end = 0ul;
  double score = 0.0;
  char strand = '\0';

  if (s >> chrom >> start >> end >> name >> score >> strand)
    region = GenomicRegion(chrom, start, end, name, score, strand);
  else region = GenomicRegion();

  char c;
  while ((c = s.get()) != '\n' && s);

  if (c != '\n')
    s.setstate(std::ios::badbit);

  if (s.eof())
    s.setstate(std::ios::badbit);

  return s;
}

bool
GenomicRegion::contains(const GenomicRegion& other) const {
  return chrom == other.chrom && start <= other.start && other.end <= end;
}


bool
GenomicRegion::overlaps(const GenomicRegion& other) const {
  return chrom == other.chrom &&
    ((start < other.end && other.end <= end) ||
     (start <= other.start && other.start < end) ||
     other.contains(*this));
}


size_t
GenomicRegion::distance(const GenomicRegion& other) const {
  if (chrom != other.chrom)
    return std::numeric_limits<size_t>::max();
  else if (overlaps(other) || other.overlaps(*this))
    return 0;
  else return (end < other.start) ?
         other.start - end + 1 : start - other.end + 1;
}


bool
GenomicRegion::operator<(const GenomicRegion& rhs) const {
  return ((chrom == rhs.chrom &&
           (start < rhs.start ||
            (start == rhs.start &&
             (end < rhs.end ||
              (end == rhs.end &&
               (strand < rhs.strand
                // || (strand == rhs.strand && name < rhs.name)
                )))))) ||
          get_chrom() < rhs.get_chrom());
}


bool
GenomicRegion::less1(const GenomicRegion& rhs) const {
  return ((chrom == rhs.chrom &&
           (end < rhs.end ||
            (end == rhs.end &&
             (start < rhs.start ||
              (start == rhs.start &&
               (strand < rhs.strand
                // || (strand == rhs.strand && name < rhs.name)
                )))))) ||
          get_chrom() < rhs.get_chrom());
}


bool
GenomicRegion::operator<=(const GenomicRegion& rhs) const {
  return !(rhs < *this);
}


bool
GenomicRegion::operator==(const GenomicRegion& rhs) const {
  return (chrom == rhs.chrom && start == rhs.start && end == rhs.end &&
          name == rhs.name && score == rhs.score && strand == rhs.strand);
}


bool
GenomicRegion::operator!=(const GenomicRegion& rhs) const {
  return (chrom != rhs.chrom || start != rhs.start || end != rhs.end ||
          name != rhs.name || score != rhs.score || strand != rhs.strand);
}


void
separate_chromosomes(const vector<SimpleGenomicRegion>& regions,
                     vector<vector<SimpleGenomicRegion> >& separated_by_chrom) {
  typedef unordered_map<chrom_id_type, vector<SimpleGenomicRegion> > Separator;
  Separator separator;
  for (vector<SimpleGenomicRegion>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const chrom_id_type the_chrom(i->chrom);
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<SimpleGenomicRegion>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}


void
separate_chromosomes(const vector<GenomicRegion>& regions,
                     vector<vector<GenomicRegion> >& separated_by_chrom) {
  typedef unordered_map<chrom_id_type, vector<GenomicRegion> > Separator;
  Separator separator;
  for (vector<GenomicRegion>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const chrom_id_type the_chrom(i->chrom);
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<GenomicRegion>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}




static bool
is_header_line(const string& line) {
  static const char *browser_label = "browser";
  static const size_t browser_label_len = 7;
  for (size_t i = 0; i < browser_label_len; ++i)
    if (line[i] != browser_label[i])
      return false;
  return true;
}


static bool
is_track_line(const char *line) {
  static const char *track_label = "track";
  static const size_t track_label_len = 5;
  for (size_t i = 0; i < track_label_len; ++i)
    if (line[i] != track_label[i])
      return false;
  return true;
}


void
ReadBEDFile(string filename, vector<GenomicRegion> &the_regions) {
  static const size_t buffer_size = 10000; // Magic

  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + filename);
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    if (!is_header_line(buffer) && !is_track_line(buffer)) {
      // const string line(buffer);
      the_regions.push_back(GenomicRegion(buffer));
    }
    in.peek();
  }
  in.close();
}


void
ReadBEDFile(string filename, vector<SimpleGenomicRegion> &the_regions) {
  std::ifstream in(filename.c_str());
  if (!in.good())
    throw BEDFileException("cannot open input file " + filename);
  size_t begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  size_t end_pos = in.tellg();
  in.seekg(0, std::ios_base::beg);

  size_t filesize = end_pos - begin_pos;
  char *buffer = new char[filesize + 1];

  in.read(buffer, filesize);
  in.close();

  the_regions.reserve(std::count(buffer, buffer + filesize, '\n'));

  char *buffer_end = buffer + filesize;
  char *c = buffer;
  while (c != buffer_end) {
    char *next = std::find(c, buffer_end, '\n');
    if (next == c)
      return;
    if (*c != 't')
      the_regions.push_back(SimpleGenomicRegion(c, next - c));
    c = next + 1;
  }
  delete[] buffer;
}


static const char *digits = "0987654321";
static const char *whitespace = " \t";

// define what to do if parts are missing: if no ':' or '-', is whole
// thing chrom?

void
parse_region_name(string region_name,
                  string& chrom, size_t &start, size_t &end) {

  const size_t colon_offset = region_name.find(":");

  // get the chromosome
  size_t chr_offset = region_name.find_last_of(whitespace, colon_offset);
  if (chr_offset == string::npos)
    chr_offset = 0;
  else
    chr_offset += 1;
  chrom = region_name.substr(chr_offset, colon_offset - chr_offset);

  // get the start
  const size_t start_end = region_name.find("-", colon_offset + 1);
  const string start_string = region_name.substr(colon_offset + 1,
                                                 start_end - colon_offset + 1);
  start = static_cast<size_t>(atoi(start_string.c_str()));

  // get the end
  const size_t end_end =
    region_name.find_first_not_of(digits, start_end + 1);
  const string end_string = region_name.substr(start_end + 1,
                                               end_end - start_end - 1);
  end = static_cast<size_t>(atoi(end_string.c_str()));
}

static size_t adjust_start_pos(const size_t orig_start,
                               const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  const size_t name_offset = chrom_name.length() + 2; // For the '>' and the '\n';
  const size_t preceding_newlines =
    static_cast<size_t>(std::floor(orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

static size_t
adjust_region_size(const size_t orig_start,
                   const string &chrom_name, const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start =
    static_cast<size_t>(std::floor(orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end =
    static_cast<size_t>(std::floor((orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}


void
extract_regions_chrom_fasta(const string &chrom_name,
                            const string &filename,
                            const vector<SimpleGenomicRegion> &regions,
                            vector<string> &sequences) {

  std::ifstream in(filename.c_str());
  for (vector<SimpleGenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size =
      adjust_region_size(orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(buffer, buffer + region_size,
                   std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    std::transform(sequences.back().begin(), sequences.back().end(),
                   sequences.back().begin(), std::ptr_fun(&toupper));
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

void
extract_regions_chrom_fasta(const string &chrom_name,
                            const string &filename,
                            const vector<GenomicRegion> &regions,
                            vector<string> &sequences) {

  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size =
      adjust_region_size(orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
                   buffer, buffer + region_size,
                   std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    std::transform(sequences.back().begin(), sequences.back().end(),
                   sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

void
extract_regions_fasta(const string &dirname,
                      const vector<GenomicRegion> &regions_in,
                      vector<string> &sequences) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
      chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw DNAmethException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(chrom_name, filenames[f_idx->second], 
                                regions[i], sequences);
  }
}

void
extract_regions_fasta(const string &dirname,
                      const vector<SimpleGenomicRegion> &regions_in,
                      vector<string> &sequences) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<SimpleGenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    // GET THE RIGHT FILE
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
      chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw DNAmethException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(chrom_name, filenames[f_idx->second], 
                                regions[i], sequences);
  }
}
