#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <filesystem>
#include <random>
#include <functional> 
#include "bioparser/fasta_parser.hpp" 


using namespace std;
namespace fs = std::filesystem;

 class Sequence  //prema  uputama za korištenje bioparser modula https://github.com/rvaser/bioparser
  { 
  public:
    Sequence(const char *name, std::uint32_t name_len, const char *data,
             std::uint32_t data_len)
    {
      name_.assign(name, name_len);
      data_.assign(data, data_len);
    }

  public:
    std::string name_, data_;
  };



string get_data_directory() { //za dohvaćanje puta do direktorija data
    filesystem::path current_path = filesystem::current_path();
    filesystem::path data_directory = current_path.parent_path().parent_path() / "data";
    return data_directory.string();
}

uint64_t sequence_length(string name) {  //fja za dohvaćanje duljine očitanja iz metapodataka fasta datoteke
    size_t last_slash = name.rfind('/');
    size_t last_underscore = name.rfind('_');
    if (last_slash != string::npos && last_underscore != string::npos) {
        string start_str = name.substr(last_slash + 1, last_underscore - last_slash - 1);
        string end_str = name.substr(last_underscore + 1);
        uint64_t start = 0;
        uint64_t end = 0;
        try {
            start = stoull(start_str);
            end = stoull(end_str);
        } catch (const std::invalid_argument& ia) {
            printf("Neispravan format u fasta datoteci: %s\n", ia.what());
            return 0;
        } catch (const std::out_of_range& oor) {
            printf("Predugačak broj u fasta datoteci: %s\n", oor.what());
            return 0;
        }
        return end - start;
    }
    printf("Neispravan format u fasta datoteci!");
    return 0;
}



vector<double> get_distribution_vector(const unique_ptr<Sequence>& fragment, int k, bool referent) {
    string name = fragment->name_;
    uint64_t length;
    if(referent==false){
        length = sequence_length(name);
        length = (int) length;       
    }
    else{
        length = fragment->data_.size();
    }
    length = length - k + 1;
    map<string, double> kmer_counts;
    string sequence = fragment->data_;
    vector<double> distribution_vector;

    // inicijalizacija mape
    for (int i = 0; i < pow(4, k); i++) {
        string kmer = "";
        int x = i;
        for (int j = 0; j < k; j++) {
            char base = "ACGT"[x % 4];
            kmer = base + kmer;
            x /= 4;
        }
        kmer_counts[kmer] = 0;
    }

    // brojanje k-mera
    for (uint64_t i = 0; i < length; i++) {
        string kmer = sequence.substr(i, k);
        kmer_counts[kmer]++;
    }
    
    // racunanje distribucije
    for (auto& pair : kmer_counts) {
        pair.second /= length;
    }

    // dodavanje komponenti u vektor odgovarajućim redoslijedom
    for (int i = 0; i < pow(4, k); i++) {
        string kmer = "";
        int x = i;
        for (int j = 0; j < k; j++) {
            char base = "ACGT"[x % 4];
            kmer = base + kmer;
            x /= 4;
        }

        distribution_vector.push_back(kmer_counts[kmer]);
    }

    return distribution_vector;
}



void write_csv(const vector<reference_wrapper<vector<double>>>& references_to_distribution_vectors, const vector<string>& names, const string& filename) {
    if (references_to_distribution_vectors.size() != names.size()) {
        printf("Velicina vektora referenci je %d,a velicina vektora imena je %d\n",(int)references_to_distribution_vectors.size(),(int)names.size());
        throw runtime_error("Veličina vektora referenci i vektora imena se ne podudara.");
    }
    
    ofstream file(filename);
    if (!file.is_open()) {
    throw runtime_error("Ne mogu otvoriti datoteku: " + filename);
}
    
    // ispis zaglavlja
    for (size_t i = 0; i < names.size(); i++) {
        file << names[i];
        if (i != names.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    size_t length = references_to_distribution_vectors[0].get().size();//to je broj redaka

    for (size_t i = 0; i < length; i++) { //redak
        for (size_t j = 0; j < references_to_distribution_vectors.size(); j++) { //stupac           
                file << references_to_distribution_vectors[j].get()[i];         
            if (j != references_to_distribution_vectors.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    
    file.close();
}


double get_cosinus_similarity(std::map<std::string, double> vector1, std::map<std::string, double> vector2) {
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;

    for (const auto& pair : vector1) {
        const auto& key = pair.first;
        const auto& value1 = pair.second;

        double value2 = 0.0;
        if (vector2.count(key) > 0) {
            value2 = vector2[key];
        }

        dot_product += value1 * value2;
        norm1 += value1 * value1;
    }

    for (const auto& pair : vector2) {
        const auto& value2 = pair.second;
        norm2 += value2 * value2;
    }

    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);

    if (norm1 == 0.0 || norm2 == 0.0) {
        throw invalid_argument("Norma jednog od vektora je 0, kosinusna slicnost nije definirana.");
    }

    return dot_product / (norm1 * norm2);
}

// Funkcija za izračunavanje srednje vrijednosti
double calculateMean(const vector<int>& data, int start, int count) {
    double sum = 0.0;
    for (int i = start; i < start + count; ++i) {
        sum += data[i];
    }
    return sum / count;
}

// Funkcija za izračunavanje standardne devijacije
double calculateStdDev(const vector<int>& data, int start, int count, double mean) {
    double sum = 0.0;
    for (int i = start; i < start + count; ++i) {
        sum += pow(data[i] - mean, 2);
    }
    return sqrt(sum / count);
}

// Funkcija za pretvaranje fastq datoteke u fasta pomoću awk
void convert_fastq_to_fasta(const string& fastq_file, const string& data_directory) {
    string fasta_file = regex_replace(fastq_file, regex("\\.fastq$"), ".fasta");
    string command = "awk '{if(NR%4==1) {printf(\">%s\\n\", substr($0, 2));} else if(NR%4==2) print;}' " + data_directory +"/Reads/"+ fastq_file + " > " + data_directory +"/Reads/"+ fasta_file;
    system(command.c_str());
}


// Funkcija za parsiranje imena bakterija iz naziva datoteka
std::string parse_bacteria_name(const std::string& filename) {
    // Provjeravamo sadrži li naziv datoteke očekivani format
    size_t fasta_pos = filename.find(".fasta");
    if (fasta_pos != std::string::npos) {
        // Pronalazimo prvu poziciju podvlake
        size_t first_underscore_pos = filename.find('_');
        if (first_underscore_pos != std::string::npos) {
            // Provjeravamo postoji li "reference" u nazivu
            size_t ref_pos = filename.find("_reference");
            size_t second_underscore_pos = (ref_pos != std::string::npos) ? ref_pos : fasta_pos;
            
            // Formiramo rezultat s prvim slovom prve riječi i cijelom drugom riječi
            std::string result;
            result += std::tolower(filename[0]);
            result += ' ';
            result += filename.substr(first_underscore_pos + 1, second_underscore_pos - first_underscore_pos - 1);
            return result;
        }
    }
    // Ako format nije odgovarajući, vraćamo neizmijenjeni naziv datoteke
    return filename;
}


// Funkcija za kreiranje metagenomskog uzorka
void create_metagenomic_sample(const vector<string>& read_files, 
                               const vector<int>& read_counts, 
                               const string& output_file,const string& data_directory) {
    ofstream outfile(output_file);
    if (!outfile.is_open()) {
    throw runtime_error("Ne mogu otvoriti datoteku: " + output_file);
}
    for (size_t i = 0; i < read_files.size(); ++i) {      
         std::cout << read_files[i]<< endl;
    }
    for (size_t i = 0; i < read_files.size(); ++i) { 
 
        string input_file = read_files[i];
        int count = read_counts[i];
        if (input_file.find(".fastq") != string::npos) {
            convert_fastq_to_fasta(input_file, data_directory);
            input_file = regex_replace(input_file, regex("\\.fastq$"), ".fasta");
        }
        string command = "awk 'BEGIN {RS=\">\"; ORS=\"\"} NR>1 {print \">\"$0}' " + data_directory + "/Reads/" + input_file + " | head -n " + to_string(count * 2) + " >> " + output_file;
        system(command.c_str());
    }
    outfile.close();
}


int main(){

    //----parsiranje referentnih genoma------------------

    // Definiranje putanje do CONFIGURATION_file.txt
    filesystem::path current_path = filesystem::current_path();
    filesystem::path path = current_path.parent_path().parent_path();
    string config_file = path.string() + "/CONFIGURATION_file.txt";
    string data_directory = get_data_directory();

    
    // Definiranje potrebnih varijabli
    vector<string> reference_genomes;
    vector<bool> circular_genomes;
    vector<string> read_files;
    vector<int> read_counts;
    string line;
    int k;
    int frag_size;
    int overlap;

    // Učitavanje konfiguracijske datoteke
    bool read_references = false, read_reads = false, k_mer = false, frag_size_flag = false, overlap_flag = false;
    ifstream file(config_file);
    if (file.is_open()) {
        while (getline(file, line)) {
            // Preskoči prazne redove i one koji počinju s '#'
            if (line.empty() || line[0] == '#') {
                // Provjeri je li linija komentar za referentne genome, očitanja, k-mer, frag_size, ili overlap
                if (line.find("NAZIVE DATOTEKA KOJE SADRŽE REFERENTNE GENOME") != std::string::npos) {
                    read_references = true;
                    read_reads = false;
                    k_mer = false;
                    frag_size_flag = false;
                    overlap_flag = false;
                } else if (line.find("NAZIVE DATOTEKA KOJE SADRŽE OČITANJA") != std::string::npos) {
                    read_reads = true;
                    read_references = false;
                    k_mer = false;
                    frag_size_flag = false;
                    overlap_flag = false;
                } else if (line.find("UNESITE ŽELJENU DULJINU K-MERA") != std::string::npos) {
                    read_reads = false;
                    read_references = false;
                    k_mer = true;
                    frag_size_flag = false;
                    overlap_flag = false;
                } else if (line.find("UNESITE ŽELJENU DULJINU ODSJEČKA REFERENTNOG GENOMA") != std::string::npos) {
                    read_reads = false;
                    read_references = false;
                    k_mer = false;
                    frag_size_flag = true;
                    overlap_flag = false;
                } else if (line.find("UNESITE DULJINU PREKLAPANJA IZMEĐU PRETHODNO ODREĐENIH ODSJEČAKA") != std::string::npos) {
                    read_reads = false;
                    read_references = false;
                    k_mer = false;
                    frag_size_flag = false;
                    overlap_flag = true;
                }
                continue;
            }

            try {
                if (read_references) {
                    // Učitavanje referentnih genoma
                    istringstream iss(line);
                    string file_name, circular_str;
                    if (!(iss >> file_name >> circular_str)) {
                        throw runtime_error("Invalid format for reference genomes");
                    }
                    bool circular = (circular_str == "true");
                    reference_genomes.push_back(file_name);
                    circular_genomes.push_back(circular);
                } else if (read_reads) {
                    // Učitavanje očitanja i broja očitanja
                    istringstream iss(line);
                    string file_name;
                    int count;
                    if (!(iss >> file_name >> count)) {
                        throw runtime_error("Invalid format for read files");
                    }
                    read_files.push_back(file_name);
                    read_counts.push_back(count);
                } else if (k_mer) {
                    // Ucitavanje željene duljine k-mera
                    istringstream iss(line);
                    if (!(iss >> k)) {
                        throw runtime_error("Invalid format for k-mer length");
                    }
                } else if (frag_size_flag) {
                    // Ucitavanje željene duljine odsjeka referentnog genoma
                    istringstream iss(line);
                    if (!(iss >> frag_size)) {
                        throw runtime_error("Invalid format for fragment size");
                    }
                } else if (overlap_flag) {
                    // Ucitavanje duljine preklapanja izmedu odsjecaka
                    istringstream iss(line);
                    if (!(iss >> overlap)) {
                        throw runtime_error("Invalid format for overlap length");
                    }
                }
            } catch (const exception &e) {
                cerr << "Error reading line: " << line << ". " << e.what() << endl;
            }
        }
        file.close();
    } else {
        cerr << "Ne mogu otvoriti konfiguracijsku datoteku." << endl;
    }

    // Ispis učitanih podataka za provjeru
    cout << "Referentni genomi:" << endl;
    for (size_t i = 0; i < reference_genomes.size(); ++i) {
        std::cout << reference_genomes[i] <<  endl;
    }
    cout << endl;
    cout << "Očitanja i broj očitanja:" << endl;
    for (size_t i = 0; i < read_files.size(); ++i) {
        std::cout << read_files[i] << " " << read_counts[i] << endl;
    }
    cout << "Željena duljina k-mera: " << k << endl;
    cout << "Željena duljina odsjeka referentnog genoma: " << frag_size << endl;
    cout << "Duljina preklapanja: " << overlap << endl;

    // Parsiranje imena bakterija za referentne genome
    vector<string> names;
    for (const auto& genome_file : reference_genomes) {
        names.push_back(parse_bacteria_name(genome_file));
    }

    // Kreiranje metagenomskog uzorka ako datoteka ne postoji
    cout << "Metagenoski uzorak stvoren je iz datoteka:"<< endl;
    string metagenomic_sample_file = data_directory + "/" + "metagenomic_sample.fasta";
    if (!fs::exists(metagenomic_sample_file)) {
        for (size_t i = 0; i < read_files.size(); ++i) {      
         cout << read_files[i]<< endl;
    }
        create_metagenomic_sample(read_files, read_counts, metagenomic_sample_file, data_directory);
    } else {
        cout << "Datoteka metagenomic_sample.fasta već postoji. Preskačem stvaranje nove datoteke." << endl;
    }

    // Parsiranje imena bakterija za metagenomski uzorak i broja očitanja
    vector<string> names_metagenomic_sample;
    for (const auto& read_file : read_files) {
        names_metagenomic_sample.push_back(parse_bacteria_name(read_file));
    }


    // Parsiranje referentnih genoma i stvaranje distribucijskih vektora
    vector<vector<unique_ptr<Sequence>>> refs;
    for (const auto& genome_file : reference_genomes) {
        auto ref_genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(data_directory + "/Referent_genomes/"+ genome_file);
        auto ref = ref_genome_parser->Parse(-1);
        refs.push_back(move(ref));
    }

    vector<vector<double>> ref_distribution_vectors(refs.size());
    vector<reference_wrapper<vector<double>>> references_to_distribution_vectors;

    for (size_t i = 0; i < refs.size(); ++i) {
        ref_distribution_vectors[i] = get_distribution_vector(ref(refs[i][0]), k, true);
        references_to_distribution_vectors.push_back(ref_distribution_vectors[i]);
    }

    // Ispis distribucijskih vektora u CSV datoteku
    string filename = "reference_vectors.csv";
    string absolute_path = data_directory + "/" + filename;
    write_csv(references_to_distribution_vectors, names, absolute_path);

    // Oslobađanje memorije nakon ispisivanja u CSV datoteku
    ref_distribution_vectors.clear();  // Briše sve vektore iz 'ref_distribution_vectors'
    ref_distribution_vectors.shrink_to_fit();  // Smanjuje kapacitet vektora na 0

    references_to_distribution_vectors.clear();  // Briše sve reference


    //----opis referentnog genoma sa više vektora---------------------------


    vector<reference_wrapper<vector<double>>> partial_references_to_referent_distribution_vectors;
    vector<shared_ptr<vector<double>>> real_vectors3;
    vector<string> names_partial_refs;  //header csv datoteke

    for (size_t ref_idx = 0; ref_idx < refs.size(); ++ref_idx) {
        auto data = refs[ref_idx][0]->data_;
        bool is_circular = circular_genomes[ref_idx];
        if(is_circular){
            for (uint32_t start_idx = 0; start_idx < data.size(); start_idx += (frag_size - overlap)) {
                string fragment;
                if (start_idx + frag_size > data.size()) {
                    // Ako je kraj fragmenta izvan genoma, produži na početak
                    fragment = data.substr(start_idx) + data.substr(0, (start_idx + frag_size) % data.size());
                } else {
                    fragment = data.substr(start_idx, frag_size);
                }

                // Simulate creating a Sequence object from the fragment
                auto seq = make_shared<Sequence>("frag", 4, fragment.c_str(), fragment.size());

                // Convert shared_ptr to unique_ptr
                auto seq_unique = make_unique<Sequence>(*seq);

                // Stvaranje shared_ptr novog vektora za svaki fragment
                shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(seq_unique, k, true));

                //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
                real_vectors3.push_back(fragment_distribution_vector);
                // Pohrana referenci na nove vektore
                partial_references_to_referent_distribution_vectors.push_back(ref(*fragment_distribution_vector));

                // Dodavanje imena bakterije u vektor imena
                names_partial_refs.push_back(names[ref_idx] + to_string(start_idx / overlap + 1));
            }
        } else {
            for (uint32_t start_idx = 0; start_idx < data.size(); start_idx += (frag_size - overlap)) {
                string fragment;
                if (start_idx + frag_size > data.size()) {
                    // Ako je kraj fragmenta izvan genoma, uzmi samo ostatak do kraja
                    fragment = data.substr(start_idx);
                } else {
                    fragment = data.substr(start_idx, frag_size);
                }

                // Simulate creating a Sequence object from the fragment
                auto seq = make_shared<Sequence>("frag", 4, fragment.c_str(), fragment.size());

                // Convert shared_ptr to unique_ptr
                auto seq_unique = make_unique<Sequence>(*seq);

                // Stvaranje shared_ptr novog vektora za svaki fragment
                shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(seq_unique, k, true));

                //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
                real_vectors3.push_back(fragment_distribution_vector);
                // Pohrana referenci na nove vektore
                partial_references_to_referent_distribution_vectors.push_back(ref(*fragment_distribution_vector));

                // Dodavanje imena bakterije u vektor imena
                names_partial_refs.push_back(names[ref_idx] + to_string(start_idx / overlap + 1));
            }
        }
    }

    //ispis u csv datoteku
    string filename4 = "PARTIAL_reference_vectors.csv";
    string absolute_path4 = get_data_directory() + "/" + filename4;
    write_csv(partial_references_to_referent_distribution_vectors, names_partial_refs, absolute_path4);

    // Nakon što su podaci iz 'refs' vektora obrađeni i više nisu potrebni
    for (auto& ref_vector : refs) {
        for (auto& seq : ref_vector) {
            seq.reset();  // Oslobađa memoriju koju je unique_ptr držao
        }
        ref_vector.clear();  // Briše sve unique_ptr iz vektora
    }
    refs.clear();  // Briše sve vektore iz vektora vektora

    // Nakon što su podaci iz 'real_vectors3' vektora obrađeni i više nisu potrebni
    for (auto& vec : real_vectors3) {
        vec.reset();  // Smanjuje brojač reference i oslobađa memoriju ako je brojač 0
    }
    real_vectors3.clear();  // Briše sve shared_ptr iz vektora


    // Parsiranje referentnih genoma koristenih u metagenomskom uzorku
    for (std::string& file : read_files) { //samo dodajem _reference da uzme referenti genom
        size_t pos = file.find(".fasta");
        if (pos != std::string::npos) {
            file.insert(pos, "_reference");
        }
    }
    vector<vector<unique_ptr<Sequence>>> refs_from_sample;
    for (const auto& genome_file : read_files) {
        auto ref_genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(data_directory + "/Referent_genomes_for_NO_ERROR_reads/"+ genome_file);
        auto ref = ref_genome_parser->Parse(-1);
        refs_from_sample.push_back(move(ref));
    }


    //----parsiranje kreiranje distribucijskih vektora očitanja---------------------------


    auto fragment_parser =
                bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                    data_directory + "/" + "metagenomic_sample.fasta");
    auto fragments = fragment_parser->Parse(-1); //fragments je vektor unique_ptr<Sequence> objekata

    // Stvaranje vektora referenci
    vector<reference_wrapper<vector<double>>> references_to_fragments_distribution_vectors;
    vector<shared_ptr<vector<double>>> real_vectors; 
    vector<int> frag_sizes; //pohranjuje duljinu svakog očitanja u uzorku

    for(int i=0;i<fragments.size();i++){
        if(fragments[i]->data_.size()<k){ //preskačem prekratka očitanja (manja od zadanog k)
            printf("Na poziciji %d u metagenomskom uzorku je prekratko ocitanje\n",i*2+2);
            int zbroj=0;
            for(int k=0;k<read_counts.size();k++){
                zbroj += read_counts[k]; 
                if(i<=zbroj-1){
                    read_counts[k]=read_counts[k]-1; //zbog izbacivanja smanjim broj
                    break;
                }
            }
        }
        else{
        frag_sizes.push_back(fragments[i]->data_.size()); //zapišem duljinu svakog očitanja
        // Stvaranje shared_ptr novog vektora za svaki fragment
        shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(fragments[i], k, false));

        //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
        real_vectors.push_back(fragment_distribution_vector);
        // Pohrana referenci na nove vektore
        references_to_fragments_distribution_vectors.push_back(ref(*fragment_distribution_vector));
        }
    }

    //header csv datoteke
    vector<string> names_frags;
    // Popunjavanje vektora names_frags
    for (size_t i = 0; i < read_counts.size(); ++i) {
        for (int j = 0; j < read_counts[i]; ++j) {
            names_frags.push_back(names_metagenomic_sample[i]);
        }
    }

    //ispis u csv datoteku
    string filename2 = "fragments_vectors.csv";
    string absolute_path2 = get_data_directory() + "/" + filename2;
    write_csv(references_to_fragments_distribution_vectors, names_frags, absolute_path2);

    // Oslobađanje memorije
    for (auto& ref_vector : references_to_fragments_distribution_vectors) {
        ref_vector.get().clear();  // Briše sadržaj svakog vektora
    }
    references_to_fragments_distribution_vectors.clear();  // Briše sve reference iz 'references_to_fragments_distribution_vectors'

    for (auto& vec : real_vectors) {
        vec.reset();  // Smanjuje brojač reference i oslobađa memoriju ako je brojač 0
    }
    real_vectors.clear();  // Briše sve shared_ptr iz vektora


    for (auto& fragment : fragments) {
        fragment.reset();  // Oslobađa memoriju koju je unique_ptr držao
    }
    fragments.clear();  // Briše sve unique_ptr iz vektora


    //----kreiranje očitanja bez greške---------------------------


    vector<reference_wrapper<vector<double>>> references_to_NO_ERROR_fragments_distribution_vectors;
    vector<shared_ptr<vector<double>>> real_vectors2;

    int size_index = 0;
    random_device rd;
    mt19937 gen(rd());

    for (size_t ref_idx = 0; ref_idx < refs_from_sample.size(); ++ref_idx) {
        for (int i = 0; i < read_counts[ref_idx]; ++i) {
            int frag_size = frag_sizes[size_index++];
            int max_start_idx = refs_from_sample[ref_idx][0]->data_.size() - frag_size;
            uniform_int_distribution<> dis(0, max_start_idx); //iz jednolike razdiobe odabirem brojeve
            int start_idx = dis(gen);

            string fragment = refs_from_sample[ref_idx][0]->data_.substr(start_idx, frag_size);

            // simuliram kreiranje objekta Sequence od fragmenta
            auto seq = make_shared<Sequence>("frag", 4, fragment.c_str(), fragment.size());

            // konvertiranje shared_ptr u unique_ptr
            auto seq_unique = make_unique<Sequence>(*seq);

            // Stvaranje shared_ptr novog vektora za svaki fragment
            shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(seq_unique, k, true));

            //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
            real_vectors2.push_back(fragment_distribution_vector);
            // Pohrana referenci na nove vektore
            references_to_NO_ERROR_fragments_distribution_vectors.push_back(ref(*fragment_distribution_vector));
        }
    }

    //ispis u csv datoteku
    string filename3 = "NO_ERROR_fragments_vectors.csv";
    string absolute_path3 = get_data_directory() + "/" + filename3;
    write_csv(references_to_NO_ERROR_fragments_distribution_vectors, names_frags, absolute_path3);

    // Oslobađanje memorije nakon ispisivanja u CSV datoteku
    for (auto& vec : real_vectors2) {
        vec.reset();  // Smanjuje brojač reference i oslobađa memoriju ako je brojač 0
    }
    real_vectors2.clear();  // Briše sve shared_ptr iz vektora

    references_to_NO_ERROR_fragments_distribution_vectors.clear();  // Briše sve reference iz vektora



    //----ispis osnovnih statistika---------------------------



    // Pretpostavljamo da frag_sizes ima točan broj elemenata
    if (frag_sizes.size() != accumulate(read_counts.begin(), read_counts.end(), 0)) {
        cerr << "Vektor frag_sizes nema očekivani broj elemenata!" << endl;
        return 1;
    }

    int start = 0;
    cout << endl;
    // Ispis zaglavlja tablice
    cout << "                  ";
    for (const auto& name : names_metagenomic_sample) {
        cout << setw(15) << name;
    }
    cout << endl;

    // Ispis srednjih vrijednosti
    cout << "average_length    ";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        double mean = calculateMean(frag_sizes, start, read_counts[i]);
        cout << setw(15) << mean;
        start += read_counts[i];
    }
    cout << endl;

    // Resetiramo početak
    start = 0;

    // Ispis standardnih devijacija
    cout << "standard_deviation";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        double mean = calculateMean(frag_sizes, start, read_counts[i]);
        double stddev = calculateStdDev(frag_sizes, start, read_counts[i], mean);
        cout << setw(15) << stddev;
        start += read_counts[i];
    }
    cout << endl;
    // Resetiramo početak
    start = 0;

    // Ispis broja ocitanja po bakteriji
    cout << "num_of_reads      ";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        cout << setw(15) << read_counts[i];
        start += read_counts[i];
    }
    cout << endl;
    cout << endl;


    return 0;
}