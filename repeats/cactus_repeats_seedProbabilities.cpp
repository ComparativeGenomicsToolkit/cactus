#include <getopt.h>
#include <string.h>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <regex>
#include <iomanip>
#include <assert.h>

#include "junctionapi/junctionapi.h"
#include "streamfastaparser.h"
#include "common.h"

#define CONFIDENCE 0.9

/*
n: number of nodes
threshold: required certainty of collapsing the graph
*/
double edgeProbabilityToCollapseGraph(int n, double threshold) {
    return 1.0 - pow((1.0 - threshold)/(double)n, 1.0/((double)(n - 1)));
}


std::string kmerToSeed(std::string kmer, std::string seedPattern) {
    std::string seed;
    assert (kmer.length() == seedPattern.length());
    for (int i = 0; i < kmer.length(); i++) {
        if (seedPattern[i] == '0') {
            seed += 'x';
        }
        else {
            seed += kmer[i];
        }

    }
    return seed;
}


typedef struct Params {
    int k = 19; //kmer length
    int maxPathLength = 20; //terminate path after this length, regardless of weight
    double minPathMultiplicity = 5;
    int scoreThreshold = 5; //Score threshold to include kmer in output
    std::string seedPattern = "1110100110010101111";
} Params;

struct Edge {
    int id;
    int start;
    int end;
    int multiplicity;
    std::vector<int> incidentEdges;

};


class CompactedDeBruijn {
    private:
        std::map<int, Edge*> edges;
        const Params &params;
        std::string genome;
        std::ofstream *statsFile;

    public:
        CompactedDeBruijn(const std::string &graphFilename,
            const std::vector<std::string> &genomes, const Params &_params,
            std::ofstream *_statsFile);
        std::map<std::string, double> findHighMultiplicityPaths(const std::string &seedPattern, std::map<std::string, std::string> &seedToPacked);

};

CompactedDeBruijn::CompactedDeBruijn(const std::string &graphFilename, 
        const std::vector<std::string> &genomes, const Params &_params, 
        std::ofstream *_statsFile): params(_params), statsFile(_statsFile) {

    TwoPaCo::JunctionPosition end;
    TwoPaCo::JunctionPosition begin;
    ChrReader chrReader(genomes);
    TwoPaCo::JunctionPositionReader reader(graphFilename.c_str());


    chrReader.NextChr(genome);
    Edge *prevEdge = NULL;
    std::cerr << "Genome length = " << genome.length() << std::endl;
    if (reader.NextJunctionPosition(begin)) {
        while (reader.NextJunctionPosition(end)) {
            Segment *segment = new Segment(begin, end, genome[begin.GetPos() + params.k],
                    TwoPaCo::DnaChar::ReverseChar(genome[end.GetPos() - 1]));

            Edge *edge;
            if (!edges.count(segment->GetSegmentId())) {
                edge = new Edge();
                edges[segment->GetSegmentId()] = edge;
                edge->id = segment->GetSegmentId();
                edge->start = begin.GetPos();
                edge->end = end.GetPos();
            }
            else {
                edge = edges[segment->GetSegmentId()];
            }
            edge->multiplicity++;
            if (prevEdge) {
                prevEdge->incidentEdges.push_back(segment->GetSegmentId());
            }
            prevEdge = edge;
            delete segment;

            begin = end;
        }
    }
}


std::map<std::string, double> CompactedDeBruijn::findHighMultiplicityPaths(const std::string &seedPattern, std::map<std::string, std::string> &seedToPacked) {
    std::map<std::string, double> seedProbabilities;

    std::set<int> seen;


    Edge *edge;
    for (auto &kv : edges) {
        int edgeId = kv.first;
        edge = kv.second;

        if(seen.count(edge->id)) continue;

        //seen[edgeId] = true;

        std::string segment;

        //edges already seen on the path, to prevent cycles
        std::string path;
        std::set<int> pathEdges;
        int pathWeight = 0;
        while(edge->multiplicity > params.minPathMultiplicity) {
            if (pathEdges.count(edge->id)) break;

            std::string segmentString;
            if (edge->id > 0) {
                segment.assign(genome.substr(edge->start,
                        (edge->end - edge->start) + params.k));
            }
            else {
                segment.assign(TwoPaCo::DnaChar::ReverseCompliment(genome.substr(edge->start, (edge->end - edge->end) + params.k)));
            }
            std::cerr << "Appending " << segment << " to path" << std::endl;

            path.append(segment);
            pathEdges.insert(edge->id);
            pathWeight += edge->multiplicity;

            //Choose next segment to traverse
            int maxMultiplicity = 0;
            Edge *bestEdge = NULL;
            for (int edgeID : edge->incidentEdges) {
                Edge *candidateEdge = edges[edgeID];
                if (candidateEdge->multiplicity > maxMultiplicity) {
                    maxMultiplicity = candidateEdge->multiplicity;
                    bestEdge = candidateEdge;
                }
            }
            edge = bestEdge;
        }
        if (path.length() == 0) continue;
        std::cerr << "Found path " << path << std::endl;

        *statsFile << path.length() << " " << pathWeight << std::endl;
        //Go through the path and set the seed multiplicities
        for (int i = 0; i < path.length() - params.k; i++) {
            std::string seed = kmerToSeed(path.substr(i, params.k), seedPattern);
            std::string packed = seedToPacked[seed];
            seedProbabilities[packed] = edgeProbabilityToCollapseGraph(pathWeight, CONFIDENCE);
        }
        path.clear();
        pathEdges.clear();
    }
    return seedProbabilities;

}


int main(int argc, char **argv) {
    Params params;
    std::string graphFilename;
    std::string seedScoresFilename;
    std::string sequenceFilename;
    std::string lastzFilename;
    std::string statsFilename;
    int seedMultiplicity = 0;
    int pathMultiplicity = 0;
    struct option longopts[] = {
        {"graphFile", required_argument, 0, 'a'},
        {"sequenceFile", required_argument, 0, 'b'},
        {"lastzFile", required_argument, 0, 'd'},
        {"statsFile", required_argument, 0, 'j'},

        {"k", required_argument, 0, 'e'},
        {"minPathMultiplicity", required_argument, 0, 'f'},
        {"maxPathLength", required_argument, 0, 'g'},
        {"scoreThreshold", required_argument, 0, 'h'},
        {"seedPattern", required_argument, 0, 'i'},

        {"seedMultiplicity", no_argument, &seedMultiplicity, 1},
        {"pathMultiplicity", no_argument, &pathMultiplicity, 1},
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
                params.minPathMultiplicity = atof(optarg);
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
    std::ofstream *statsFile = NULL;
    if (statsFilename.length() > 0) {
        statsFile = new std::ofstream(statsFilename);
    }

    std::map<std::string, int> seedToMultiplicity;
    std::map<std::string, std::string> seedToPacked;

    //Parse the seed counts from lastz
    std::string seed;
    std::string packed;
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
        ss << match[1];
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
        seedToMultiplicity[packed] = count;
    }

    std::map<std::string, double> seedProbabilities;
    if (seedMultiplicity) {
        for (auto &kv : seedToMultiplicity) {
            std::string packed = kv.first;
            double multiplicity = (double)kv.second;
            if (multiplicity > 1.0) {
                seedProbabilities[packed] = edgeProbabilityToCollapseGraph(multiplicity, CONFIDENCE);
            }
            else {
                seedProbabilities[packed] = 1.0;
            }

        }
            
    }
    else if (pathMultiplicity) {
        std::vector<std::string> genomes;
        genomes.push_back(sequenceFilename);
        CompactedDeBruijn graph(graphFilename, genomes, params, statsFile);
        seedProbabilities = graph.findHighMultiplicityPaths(params.seedPattern, seedToPacked);
    }

	for(auto& kv : seedProbabilities) {
        std::string packed = kv.first;
        double probability = kv.second;
        printf("%s %lf\n", packed.c_str(), probability);
    }
    if (statsFile) {
        statsFile->close();
    }
}
