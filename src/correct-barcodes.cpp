#include <algorithm>
#include <cerrno>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <cmath>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;

std::map<char, char> nucleotide_complements = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'},
    {'N', 'N'},
    {'a', 't'},
    {'t', 'a'},
    {'c', 'g'},
    {'g', 'c'},
    {'n', 'n'},
};

std::string complement(const std::string& seq) {
    std::string complemented;
    for (const auto it : seq) {
        complemented += nucleotide_complements.at(it);
    }
    return complemented;
}

std::string reverse_complement(const std::string& seq) {
    std::string rc;
    for (int i = seq.length() - 1; i >= 0; i--) {
        rc += nucleotide_complements.at(seq[i]);
    }
    return rc;
}

class FileException: public std::runtime_error {
public:
    FileException() : std::runtime_error("") { }
    explicit FileException(std::string msg) : std::runtime_error(msg) { }
};

bool is_gzipped_file(const std::string& filename) {
    bool gzipped = false;
    FILE* f = fopen(filename.c_str(), "rb");

    if (f == NULL) {
        throw FileException("Could not open file \"" + filename + "\": " + strerror(errno));
    } else {
        if (fgetc(f) == 0x1f && fgetc(f) == 0x8b) {
            gzipped = true;
        }
    }
    fclose(f);
    return gzipped;
}

bool is_gzipped_filename(const std::string& filename) {
    std::string ext = ".gz";
    return std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}

void log_message(const std::string& message) {
    std::time_t now = std::time(0);
    std::tm* local = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", local);
    std::cerr << buffer << " " << message << std::endl;
}

// void log_message(const double& message) {
//     std::time_t now = std::time(0);
//     std::tm* local = std::localtime(&now);
//     char buffer[80];
//     std::strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", local);
//     std::cerr << buffer << " " << message << std::endl;
// }

class FASTQRecord {
public:
    std::string name;
    std::string comment;
    std::string sequence;
    std::string quality;

    void transform(int offset, int length, bool rc) {
        sequence = sequence.substr(offset, length);
        quality = quality.substr(offset, length);
        if (rc) {
            sequence = reverse_complement(sequence);
            std::reverse(quality.begin(), quality.end());
        }
    }

    void strip_description() {
        // the description is the part of the name after the first space
        if (name.find(' ')!=std::string::npos) {
            name = name.substr(0, name.find(' '));
        }
    }

};

std::ostream& operator<<(std::ostream& os, const FASTQRecord& record) {
    os << record.name << '\n'
       << record.sequence << '\n'
       << record.comment << '\n'
       << record.quality << '\n';
    return os;
}

std::istream& operator>>(std::istream& is, FASTQRecord& record) {
    getline(is, record.name);
    getline(is, record.sequence);
    getline(is, record.comment);
    getline(is, record.quality);
    return is;
}

class TrieNode {
public:
    std::unordered_map<char, TrieNode*> children;
    bool is_end_of_word = false;
};


class Trie {
public:
    TrieNode root;

    void insert(const std::string& word) {
        TrieNode* current = &root;
        for (const auto& letter : word) {
            if (current->children.find(letter) == current->children.end()) {
                current->children[letter] = new TrieNode();
            }
            current = current->children[letter];
        }
        current->is_end_of_word = true;
    }

    void build_trie(const std::unordered_set<std::string>& words) {
        for (const auto& word : words) {
            insert(word);
        }
    }

    int hamming_distance(const std::string& word1, const std::string& word2) {
        if (word1.size() != word2.size()) {
            throw std::runtime_error("Cannot calculate hamming distance for two words of different lengths");
        }
        return std::inner_product(word1.begin(), word1.end(), word2.begin(), 0, std::plus<int>(), std::not_equal_to<char>());
    }

    void find_similar_words(const std::string& word, int max_distance, std::vector<std::string>& result) {
        _search(&root, word, "", result, max_distance);
    }

    void _search(TrieNode* node, const std::string& word, const std::string& current_word, std::vector<std::string>& result, int max_distance, int current_distance=0) {
        if (current_distance > max_distance) {
            return;
        }
        if (current_word.size() == word.size()) {
            if (node->is_end_of_word && hamming_distance(current_word, word) <= max_distance) {
                result.push_back(current_word);
            }
            return;
        }

        for (auto& it : node->children) {
            auto& letter = it.first;
            TrieNode* child_node = it.second;
            if (current_word.size() < word.size()) {
                int next_distance = current_distance + (letter != word[current_word.size()]);
                _search(child_node, word, current_word + letter, result, max_distance, next_distance);
            } else {
                _search(child_node, word, current_word + letter, result, max_distance, current_distance + 1);
            }
        }
    }
};

std::unordered_map<std::string,long long int> read_counts(const std::string& counts_filename) {
    std::ifstream counts_file(counts_filename);
    std::unordered_map<std::string,long long int> counts;
    std::string line;

    while (std::getline(counts_file, line)) {
        // split the line
        std::istringstream iss(line);
        std::string barcode;
        long long int count;
        iss >> barcode >> count;
        counts[barcode] = count;
    }

    counts_file.close();

    return counts;
}

std::unordered_set<std::string> read_whitelist(const std::string& whitelist_filename) {
    std::ifstream whitelist_file(whitelist_filename);
    std::unordered_set<std::string> whitelist_barcodes;
    std::string line;

    while (std::getline(whitelist_file, line)) {
        whitelist_barcodes.insert(line);
    }

    whitelist_file.close();

    return whitelist_barcodes;
}

double likelihood_of_errors(const std::string& uncorrected, const std::string& potential_correction, const std::string& phred) {
    double total_p = 1.0;
    for (int i = 0; i < uncorrected.size(); i++) {
        if (uncorrected[i] != potential_correction[i]) {
            int q = phred[i] - 33;
            double p = std::pow(10, -q/10.0);
            total_p *= p;
        }
    }
    return total_p;
}

std::string correct_barcode(const std::string& uncorrected, const std::vector<std::string>& potential_fixes, const std::string& phred, std::unordered_map<std::string, long long int>& counts) {
    // log_message("Correcting barcode (phred): " + uncorrected + " (" + phred + ")");
    // log_message("Potential fixes: ");
    // for (const auto& fix : potential_fixes) {
    //     log_message(fix);
    // }
    if (potential_fixes.size() == 1) {
        return potential_fixes[0];
    } else if (potential_fixes.size() == 0) {
        return "";
    } else {
        std::vector<long long int> priors;
        for (const auto& barcode : potential_fixes) {
            priors.push_back(counts[barcode]);
        }

        std::vector<double> prob_of_errors;
        for (const auto& barcode : potential_fixes) {
            prob_of_errors.push_back(likelihood_of_errors(uncorrected, barcode, phred));
        }

        std::vector<double> relative_magnitude;
        for (int i = 0; i < priors.size(); i++) {
            relative_magnitude.push_back(priors[i] * prob_of_errors[i]);
        }

        // log_message("Prob of errors: ");
        // for (const auto& pe : prob_of_errors) {
        //     log_message(pe);
        // }
        // log_message("Priors: ");
        // for (const auto& p : priors) {
        //     log_message(std::to_string(p));
        // }
        // log_message("Relative magnitude: ");
        // for (const auto& rm : relative_magnitude) {
        //     log_message(rm);
        // }

        double norm_factor = std::accumulate(relative_magnitude.begin(), relative_magnitude.end(), 0.0);
        if (norm_factor == 0.0) {
            // log_message("No good fix found (accumulated 0).");
            return "";
        }

        for (int i = 0; i < relative_magnitude.size(); i++) {
            relative_magnitude[i] = relative_magnitude[i] / norm_factor;
        }
        // log_message("Relative magnitudes, normalized: ");
        // for (const auto& rm : relative_magnitude) {
        //     log_message(rm);
        // }

        if (*std::max_element(relative_magnitude.begin(), relative_magnitude.end()) >= 0.975) {
            for (int i = 0; i < potential_fixes.size(); i++) {
                if (relative_magnitude[i] >= 0.975) {
                    return potential_fixes[i];
                }
            }
        } else {
            // log_message("No good fix found (max relative magnitude).");
            // log_message("Max relative magnitude:");
            // log_message(*std::max_element(relative_magnitude.begin(), relative_magnitude.end()));
            return "";
        }
    }
}




void print_usage() {
    std::cout << "Correct cell barcodes in a 10X snATAC-seq library.\n\n"

              << "Usage:\n\n" << "correct-barcodes [options] input_file barcode_whitelist output_file\n\n"
              << "where:\n"
              << "    input_file is the fastq file of barcode reads\n"
              << "    barcode_whitelist is the barcode whitelist\n"
              << "    barcode_counts is the file of barcode counts (before correction)\n"
              << "    output_file will be the fastq file of barcodes\n\n"

              << "Options may include:\n\n"

              << "-h|--help: show this usage message.\n"
              << "-v|--verbose: show more details and progress updates." << std::endl << std::endl;

}

int main(int argc, char **argv)
{
    int c, option_index = 0;
    bool verbose = 0;

    std::string input_filename;
    std::string whitelist_filename;
    std::string counts_filename;
    std::string output_fastq;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", no_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // parse the command line arguments
    while ((c = getopt_long(argc, argv, "vh?", long_options, &option_index)) != -1) {
        switch (c) {
        case 'h':
        case '?':
            print_usage();
            exit(1);
        case 'v':
            verbose = true;
            break;
        case 0:
            if (long_options[option_index].flag != 0){
                break;
            }

        default:
            print_usage();
            exit(1);
        }
    }

    if (argc != (optind + 4)) {
        print_usage();
        exit(1);
    }

    input_filename = argv[optind];
    whitelist_filename = argv[optind + 1];
    counts_filename = argv[optind + 2];
    output_fastq = argv[optind + 3];

    std::unordered_set<std::string> whitelist_barcodes = read_whitelist(whitelist_filename);

    if (verbose) {
        log_message("Read " + std::to_string(whitelist_barcodes.size()) + " barcodes from whitelist.");
    }

    if (verbose) {
        log_message("Reading counts...");
    }
    std::unordered_map<std::string, long long int> counts = read_counts(counts_filename);
    
    // filter counts to barcodes on the whitelist
    std::unordered_map<std::string, long long int> whitelist_counts;
    for (const auto& barcode : whitelist_barcodes) {
        if (counts.find(barcode) != counts.end()) {
            whitelist_counts[barcode] = counts[barcode];
        }
    }

    if (verbose) {
        log_message("Read " + std::to_string(counts.size()) + " barcodes from counts.");
        log_message("Filtered counts to " + std::to_string(whitelist_counts.size()) + " barcodes in whitelist.");
    }

    // add a pseudocount to all whitelist barcodes
    for (const auto& barcode : whitelist_barcodes) {
        if (whitelist_counts.find(barcode) == whitelist_counts.end()) {
            whitelist_counts[barcode] = 1;
        } else {
            whitelist_counts[barcode]++;
        }
    }

    if (verbose) {
        log_message("Building trie...");
    }
    Trie trie;
    trie.build_trie(whitelist_barcodes);

    if (verbose) {
        log_message("Finished building trie...");
    }

    if (verbose) {
        log_message("Reading counts...");
        std::unordered_map<std::string, long long int> counts;
    }

    // input stream
    bio::file_source input_file(input_filename);
    bio::stream<bio::file_source> input_stream(input_file);

    bio::filtering_stream<bio::input> in;
    if (is_gzipped_file(input_filename)) {
        in.push(bio::gzip_decompressor());
    }
    in.push(input_stream);


    //output fastq
    std::ofstream output_file(output_fastq, std::ios_base::out | std::ios_base::binary);

    bio::filtering_stream<bio::output> out;
    if (is_gzipped_filename(output_fastq)) {
        out.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    }
    out.push(output_file);


    FASTQRecord record;
    long long int record_count = 0;
    long long int number_already_in_whitelist = 0;
    long long int number_corrected = 0;
    long long int number_not_corrected = 0;

    std::unordered_map<std::string,std::vector<std::string>> potential_corrections;

    while (in >> record) {
        record_count++;
        if (record_count % 1000000 == 0 && verbose) {
            log_message("Processed " + std::to_string(record_count) + " records so far...");
        }

        record.strip_description();

        if (whitelist_barcodes.find(record.sequence) != whitelist_barcodes.end()) {
            record.name = record.name + " CR:Z:" + record.sequence + "\tCB:Z:" + record.sequence + "\tCY:Z:" + record.quality;
            out << record;
            number_already_in_whitelist++;
        } else {

            std::vector<std::string> similar;

            if (potential_corrections.find(record.sequence) == potential_corrections.end()) {   
                trie.find_similar_words(record.sequence, 2, similar);
                potential_corrections[record.sequence] = similar;
            } else {
                similar = potential_corrections[record.sequence];
            }

            std::string corrected = correct_barcode(record.sequence, similar, record.quality, whitelist_counts);

            if (corrected != "") {
                record.name = record.name + " CR:Z:" + record.sequence + "\tCB:Z:" + corrected + "\tCY:Z:" + record.quality;
                record.sequence = corrected;
                number_corrected++;
            } else {
                record.name = record.name + " CR:Z:" + record.sequence + "\tCY:Z:" + record.quality;
                number_not_corrected++;
            }
            out << record;
        }

    }

    bio::close(in);
    bio::close(out);

    if (verbose) {
        log_message("Processed " + std::to_string(record_count) + " records.");
        log_message("Number of records already in whitelist: " + std::to_string(number_already_in_whitelist));
        log_message("Number of records corrected: " + std::to_string(number_corrected));
        log_message("Number of records not corrected: " + std::to_string(number_not_corrected));
        log_message("Done.");
    }

}
