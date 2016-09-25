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


class CompactedVertex {
    private:
        int position;
        std::map<int, int> edges;
        int multiplicity;
    public:
        void addEdge(int toId);
        int stepForward();
        int getPosition() {return position;}
        int getMultiplicity() {return multiplicity;}
        CompactedVertex(int _position): position(_position) {}
};

int CompactedVertex::stepForward() {
    int highestMultiplicity = 0;
    int bestEdge = 0;
    for (auto &kv : edges) {
        if (kv.second > highestMultiplicity) {
            highestMultiplicity = kv.second;
            bestEdge = kv.first;
        }
    }
    return bestEdge;
}
void CompactedVertex::addEdge(int id) {
    edges[id] += 1;
    multiplicity += 1;
}

class CompactedDeBruijn {
    private:
        std::map<int, Vertex*> vertices;
        std::string sequence;
        const Params &params;

    public:
        CompactedDeBruijn(const std::string &sequenceFilename, 
                const std::string &graphFilename,
                const Params &_params);
        std::map<std::string, int> findHighMultiplicityPaths();

};

CompactedDeBruijn::CompactedDeBruijn(const std::string &sequenceFilename,
        const std::string &graphFilename): params(_params) {
    params = _params;
    TwoPaCo::JunctionPosition pos;
    TwoPaCo::JunctionPositionReader reader(graphFilename.c_str());


    reader.NextJunctionPosition(pos);
    int prevId = pos.GetId();
    vertices[prevId] = new Vertex(pos.GetPos());
    //Index junction IDs by position and position by junction ID
    while (reader.NextJunctionPosition(pos)) {
        int id = pos.GetId();
        if (!vertices.count(id)) {
            vertices[id] = new Vertex(pos.GetPos());
        }
        vertices[prevId].addEdge(id);


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

std::map<std::string, int> findHighMultiplicityPaths() {
    std::set<int> seen;
    std::vector<int> path;
    std::map<std::string, int> kmerMultiplicity;

    Vertex* v = NULL;
    for (auto &kv : vertices) {
        int id = kv.first;
        v = kv.second;
        if(seen.find(v) != v.back()) continue;

        //Try to extend a high-multiplicity path starting at
        //this vertex
        int pathWeight = 0;
        int multiplicity = v->getMultiplicity();
        while (multiplicity > params.minPathMultiplicity) {

            //Don't traverse cycles
            if(std::find(path, id) != path.back()) break;

            path.push_back(id);
            seen.insert(id);
            pathWeight += multiplicity;


            id = v->stepForward();
            v = vertices[id];
            multiplicity = v->getMultiplicity();

        }
        if (path.size() == 0) continue;

        int pathStart = vertices[path[0]]->getPosition();
        int pathStop = vertices[path.back()]->getPosition();
        for (int i = pathStart; i < pathStop + params.k; i++) {
            kmerMultiplicity[sequence.substr(i, params.k)] = pathWeight;
        }
        path.clear();

    }
}

std::string kmerToSeed(std::string kmer, std::string seedPattern) {
    std::string seed;
    assert (kmer.length() == seedPattern.length());
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
    std::map<unsigned int, int> seedMultiplciity = parseLastzTable(lastzFilename);
    SeedProbabilities seedProbabilities;
    if (seedMultiplicity) {
        seedProbabilities = seedProbabilitiesFromSeedMultiplicity(seedMultiplicity);
    }
    else if (highMultiplicityPath) {
        CompactedDeBruijn graph(sequenceFilename, graphFilename, params);
        seedProbabilities = seedProbabilitiesFromPathMultiplicity(seedMultiplicity, graph);
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
