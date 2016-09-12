#include <getopt.h>
#include <string.h>
#include <vector>
#include <map>
#include <iostream>
#include <regex>
#include <iomanip>
#include <assert.h>

#include "junctionapi/junctionapi.h"
#include "streamfastaparser.h"


typedef std::map<int, int> JunctionMultiplicity;
typedef std::map<int, int> PositionToJunction;
typedef std::map<int, std::vector<int>> JunctionPositions;
typedef std::map<std::string, double> SeedProbabilities;


/*
n: number of nodes
threshold: required certainty of collapsing the graph
*/
double edgeProbabilityToCollapseGraph(int n, double threshold) {
    return 1.0 - pow((1.0 - threshold)/(double)n, 1.0/((double)(n - 1)));
}

typedef struct Params {
    int k = 19; //kmer length
    int maxPathLength = 20; //terminate path after this length, regardless of weight
    double pathThreshold = 1.25; //threshold of pathWeight/pathLength at which to terminate path
    int scoreThreshold = 5; //Score threshold to include kmer in output
    std::string seedPattern = "1110100110010101111";
} Params;


class SeedSampler {
    private:
        Params params;
        std::string sequence;

        JunctionPositions junctionPositions;
        PositionToJunction positionToJunction;
        std::ofstream *statsFile = NULL;

    public:
        int getJunctionID(int position);
        int extendPath(bool forward, std::vector<int> *positions);
        int getMultiplicity(int position);
        std::vector<int> *getPositions(int position);
        SeedProbabilities sampleHighMultiplicityPaths();
        std::string kmerToSeed(std::string kmer);
        SeedProbabilities seedMultiplicity();
        SeedSampler(const std::string &sequenceFilename, 
                const std::string &graphFilename, 
                const std::string &statsFilename,
                const Params &params);
        ~SeedSampler();
                


};

SeedSampler::SeedSampler(const std::string &sequenceFilename,
        const std::string &graphFilename, 
        const std::string &statsFilename, 
        const Params &_params) {

    if (!statsFilename.empty()) {
        statsFile = new std::ofstream(statsFilename);
    }

    params = _params;
    TwoPaCo::JunctionPosition pos;
    TwoPaCo::JunctionPositionReader reader(graphFilename.c_str());


    //Index junction IDs by position and position by junction ID
    while (reader.NextJunctionPosition(pos)) {
        positionToJunction[pos.GetPos()] = pos.GetId();
        junctionPositions[pos.GetId()].push_back(pos.GetPos());

    }

    //Read the genome
    TwoPaCo::StreamFastaParser parser(sequenceFilename);
    char ch;
    if (parser.ReadRecord()) {
        while (parser.GetChar(ch)) {
            sequence.push_back(ch);
        }
    }


}

SeedSampler::~SeedSampler() {
    if (statsFile) {
        statsFile->close();
        delete statsFile;
    }
}

std::vector<int> * SeedSampler::getPositions(int position) {
    if(positionToJunction.count(position)) {
        int junctionID = positionToJunction[position];
        return &(junctionPositions[junctionID]);
    }
    return NULL;

}

int SeedSampler::getMultiplicity(int position) {
    if (positionToJunction.count(position)) {
        int id = positionToJunction[position];
        return junctionPositions[id].size();
    }
    else {
        return 1;
    }
}


int SeedSampler::extendPath(bool forward, std::vector<int> *positions) {

    //Find the highest multiplicity adjacent vertex in the given direction
    int bestPosition = 0;
    int maxMultiplicity = 0;
    for (int position : *positions) {
        int adjPos = (forward) ? position + 1 : position - 1;
        int adjMultiplicity = this->getMultiplicity(adjPos);
        if (adjMultiplicity > maxMultiplicity) {
            maxMultiplicity = adjMultiplicity;
            bestPosition = adjPos;
        }
    }
    return bestPosition;
}

SeedProbabilities SeedSampler::sampleHighMultiplicityPaths() {

    //IDs of junctions that have already been included in a path
    std::vector<int> seen;

    //Integrated multiplicity of the path
    //that contains each kmer
    SeedProbabilities seedProbabilities;


    //Try to find a high multiplicity path starting at each junction
    for (auto &kv : positionToJunction) {
        int junctionID = kv.second;
        int multiplicity = getMultiplicity(junctionID);

        //list of positions in this path
        std::vector<int> path;
        std::vector<int> *positions = &(junctionPositions[junctionID]);

        //extend the path in both directions by moving to the highest
        //multiplicity adjacent position, unless the path becomes too "thin"
        double pathWeight = (double)multiplicity;
        double pathLength = 1.0;

        while((pathWeight/pathLength > params.pathThreshold) && (pathLength < params.maxPathLength)) {
            //Add the current junction to the path, using the first
            //position in the list as its representative
            path.push_back(positions->at(0));

            int newPosition = this->extendPath(true, positions);
            pathWeight += (double)getMultiplicity(newPosition);
            pathLength += 1.0;
            positions = getPositions(newPosition);
            if (!positions) {
                positions = new std::vector<int>();
                positions->push_back(newPosition);
            }
        }
        std::cerr << "Found path with length " << pathLength << " and weight "
            << pathWeight << std::endl;
        assert(statsFile);
        if(statsFile) {
            *statsFile << pathLength << "\t" << pathWeight << std::endl;
        }


        //Go through the path and set the score of each kmer
        //to the integrated multiplicity of the path
        for (int position : path) {
            std::string kmer = sequence.substr(position, params.k);
            std::string seed = this->kmerToSeed(kmer);
            if (pathWeight > params.scoreThreshold) {
                seedProbabilities[seed] = edgeProbabilityToCollapseGraph(pathWeight, 0.9);
            }

        }


    }
    return seedProbabilities;
}

SeedProbabilities SeedSampler::seedMultiplicity() {
    SeedProbabilities seedProbs;
    for (auto &kv : junctionPositions) {
        int junctionID = kv.first;
        std::vector<int> positions = kv.second;
        int multiplicity = positions.size();
        std::string kmer = sequence.substr(positions[0], params.k);
        std::string seed = this->kmerToSeed(kmer);
        seedProbs[seed] = edgeProbabilityToCollapseGraph(multiplicity, 0.9);
    }
    return seedProbs;

}

std::string SeedSampler::kmerToSeed(std::string kmer) {
    std::string seed;
    for (int i = 0; i < kmer.length(); i++) {
        if (params.seedPattern[i] == '0') {
            seed += 'x';
        }
        else {
            seed += kmer[i];
        }

    }
    return seed;
}

std::map<std::string, unsigned long> parseLastzTable(const std::string &lastzFilename) {
    std::map<std::string, unsigned long> seedToPacked;
    std::string seed;
    unsigned long packed;
    unsigned int count;

    std::ifstream lastzFile(lastzFilename);
    std::string line;
    const std::regex pattern("([0-9A-F]+)\\/([ACTGx]+):\\s+([0-9]+)");
    std::stringstream ss;
    std::smatch match;
    while (std::getline(lastzFile, line)) {
        std::regex_match(line, match, pattern);
        if (match.size() < 3) continue;
        ss.clear();
        ss << std::hex << match[1];
        ss >> packed;
        ss.clear();
        ss << match[2];
        ss >> seed;
        ss.clear();
        ss << match[3];
        ss >> count;
        ss.clear();

        //std::cerr << "Parsed seed " << seed << " with packed representation " 
        //    << packed << std::endl;
        seedToPacked[seed] = packed;
    }
    return seedToPacked;

}



int main(int argc, char **argv) {
    Params params;
    std::string graphFilename;
    std::string seedScoresFilename;
    std::string sequenceFilename;
    std::string lastzFilename;
    std::string statsFilename;
    int seedMultiplicity = 0;
    int highMultiplicityPath = 0;
    struct option longopts[] = {
        {"graphFile", required_argument, 0, 'a'},
        {"sequenceFile", required_argument, 0, 'b'},
        {"lastzFile", required_argument, 0, 'd'},
        {"statsFile", required_argument, 0, 'j'},

        {"k", required_argument, 0, 'e'},
        {"pathThreshold", required_argument, 0, 'f'},
        {"maxPathLength", required_argument, 0, 'g'},
        {"scoreThreshold", required_argument, 0, 'h'},
        {"seedPattern", required_argument, 0, 'i'},

        {"seedMultiplicity", no_argument, &seedMultiplicity, 1},
        {"pathMultiplicity", no_argument, &highMultiplicityPath, 1},
        {0, 0, 0, 0} 
    };
    int flag;
    while((flag = getopt_long(argc, argv, "a:b:d:e:f:g:h:i:j:", longopts, NULL)) != -1) {
        switch(flag) {
            case 'a':
                graphFilename.assign(optarg);
                break;
            case 'b':
                sequenceFilename.assign(optarg);
                break;
            case 'd':
                lastzFilename.assign(optarg);
                break;
            case 'j':
                statsFilename.assign(optarg);
                break;
            case 'e':
                params.k = atoi(optarg);
                break;
            case 'f':
                params.pathThreshold = atof(optarg);
                break;
            case 'g':
                params.maxPathLength = atoi(optarg);
                break;
            case 'h':
                params.scoreThreshold = atoi(optarg);
                break;
            case 'i':
                params.seedPattern = std::string(optarg);
                break;
            case '?':
                break;
            default:
                break;
        }
    }
    SeedSampler seedSampler(sequenceFilename, graphFilename, statsFilename, params);
    SeedProbabilities seedProbabilities;
    if (seedMultiplicity) {
        seedProbabilities = seedSampler.seedMultiplicity();
    }
    else if (highMultiplicityPath) {
        seedProbabilities = seedSampler.sampleHighMultiplicityPaths();
    }


    //Get the packed representation of each seed from Lastz's position table for now
    std::map<std::string, unsigned long> seedToPacked = parseLastzTable(lastzFilename);


	for(auto& kv : seedToPacked) {
        std::string seed = kv.first;
        unsigned long packed = kv.second;
        double probability = 1.0;
        if (seedProbabilities.count(seed)) {
            probability = seedProbabilities[seed];
        }
        printf("%lx %f\n", packed, probability);
    }
}
