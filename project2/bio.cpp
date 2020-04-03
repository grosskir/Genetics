 #include <algorithm>
 #include <iostream>
 #include <sstream>
 #include <stdexcept>
 #include <string>
 #include <vector>
 bool IsValidDNASequence(const std::string& input) {
     std::istringstream input_istream(input);
     char c_in;
     bool valid = true;

     while (input_istream >>std::noskipws>> c_in) {

         if (c_in == 'A' || c_in == 'T' || c_in == 'C' || c_in == 'G') {
             valid = true;
         } else {
             valid = false;
             break;
         }
     }
     return valid;
 }

std::string switch_chars(char value, const std::string& input) {
     std::istringstream input_istream(input);
     std::ostringstream output_stream;
     char c_in;
     std::string output;

     while (input_istream >> c_in) {
         if (c_in == 'C') {
             output_stream << 'G';
         } else if (c_in == 'G') {
             output_stream << 'C';
         } else if (c_in == 'A') {
             output_stream << value;
         } else {
             output_stream << 'A';
         }
     }
     output = output_stream.str();
     std::reverse(output.begin(), output.end());
     return output;
 }

 void GetReverseComplementSequence(const std::string& input,
                                   std::string* const output) {
     *output = switch_chars('T', input);
 }
 std::string GetRNATranscript(const std::string& input) {
     return switch_chars('U', input);
 }

 std::vector<std::vector<std::string>> threeReadingFrames(
     const std::vector<char>& input) {
     // int index = 0;
     std::vector<std::vector<std::string>> all_sequences;
     std::vector<std::string> sequence;
     std::string codon;
     for (int shift = 0; shift < 3; shift++) {
         for (int counter = shift; counter < static_cast<int>(input.size());
              counter += 3) {
             for (int codon_val = 0; codon_val < 3; codon_val++) {
                 try {
                     codon.push_back(input.at(counter + codon_val));
                 } catch (std::out_of_range& e) {
                     break;
                 }
             }
             if (static_cast<int>(codon.size()) == 3) {
                 sequence.push_back(codon);
             }
             codon = "";
         }

         all_sequences.push_back(sequence);
         sequence = {};
         // index = shift;
     }
     return all_sequences;
 }
std::vector<std::vector<std::string>> GetReadingFramesAsCodons(
     const std::string& input) {
     std::string RNA = GetRNATranscript(input);
     std::string antiRNA = GetRNATranscript(RNA);
     // std::vector<char> original
     // =static_cast<std::vector<char>>(GetRNATranscript(input));
     // std::vector<char> antiparallel(original.end(),original.begin());
     std::vector<std::vector<std::string>> all_sequences;
     std::vector<char> original(RNA.begin(), RNA.end());
     std::vector<char> antiparallel(antiRNA.begin(), antiRNA.end());



     all_sequences = threeReadingFrames(original);
     std::vector<std::vector<std::string>> anti_codons = threeReadingFrames(antiparallel);
     for(const auto codon:anti_codons){
         all_sequences.push_back(codon);
     }

     return all_sequences;
 }

std::string Translate(const std::vector<std::string>& codon_sequence) {
     std::vector<std::string> possible_codons = {
         "GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
         "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
         "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
         "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
         "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
         "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
         "GUG", "UAG", "UGA", "UAA"};
     std::vector<std::string> amino_acids = {
         "A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D",
         "D", "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H",
         "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F",
         "F", "P", "P", "P", "P", "S", "S", "S", "S", "S", "S", "T", "T",
         "T", "T", "W", "Y", "Y", "V", "V", "V", "V", "*", "*", "*"};

     std::string translated_codons = "";
     std::string one_amino_acid;

     for (const std::string codon_val : codon_sequence) {
	 // from https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/
         auto p = std::find(std::begin(possible_codons),
                            std::end(possible_codons), codon_val);
         auto d = std::distance(std::begin(possible_codons), p);
         one_amino_acid = amino_acids[d];
         translated_codons += one_amino_acid;
     }

     return translated_codons;
 }

std::string GetLongestOpenReadingFrame(const std::string& DNA_sequence) {
     std::vector<std::vector<std::string>> all_reading_frames =
         GetReadingFramesAsCodons(DNA_sequence);
     std::string longest_reading_frame = "";
     std::string translated;
     std::ostringstream test_stream;
     char amino_acid;
     bool valid_stream = false;
     std::string test_string;
     for (const auto reading_frame : all_reading_frames) {
         translated = Translate(reading_frame);
         std::istringstream frame_istream(translated);
        
         test_stream.str(std::string());
         while (frame_istream >> amino_acid) {
             if (amino_acid == 'M') {
                 valid_stream = true;
             }
             if (valid_stream == true) {
                 test_stream << amino_acid;
             }
             if (amino_acid == '*') {
                 valid_stream = false;
                 test_string = test_stream.str();
               
                 if (static_cast<int>(longest_reading_frame.size()) <
                     static_cast<int>(test_string.size())) {
                     longest_reading_frame = test_string;
                 }
                 test_stream.str(std::string());
                 test_string = "";
             }
         }
     }
     return longest_reading_frame;
 }


